"""
parte1_base.py
Parte 1 do TC5: simulação base (linear, mur = 1200) na posição nominal x = X_GAP.

O que o script faz:
  1. Resolve o problema com I = 1 A.
  2. Salva figuras: malha, A' com linhas de fluxo, |B|, |H| e Hr/Hz na esfera.
  3. Calcula Fz na esfera de dois jeitos (Maxwell manual e FacetNormal+assemble).
  4. Calcula a corrente de equilíbrio (Fz = peso), refaz a simulação nela e confere.
  5. Verifica o fio AWG 18 (densidade de corrente e potência dissipada).
  6. Salva tudo em results/parte1/ (PNG, CSV, JSON e resumo em texto).
     O JSON é lido pela Parte 2 (corrente corrigida).

Uso:
    python parte1_base.py
    python parte1_base.py --x-gap-mm 5 --current 1
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from dolfinx import fem

import config as cfg
import geometry
import solver
import postprocess
import viz

OUT = Path(__file__).resolve().parent / "results" / "parte1"


def dados_triangulacao(sol):
    """Pontos, triângulos e A' nodal (P1) para os gráficos em matplotlib."""
    V1 = fem.functionspace(sol.mesh, ("Lagrange", 1))
    pontos = V1.tabulate_dof_coordinates()[:, :2]
    lista = V1.dofmap.list
    if callable(lista):          # compatibilidade entre versões do dolfinx
        lista = lista()
    tri = np.asarray(lista).reshape(-1, 3)
    A1 = fem.Function(V1)
    A1.interpolate(sol.A)        # P2 -> P1, só para desenhar
    return pontos, tri, A1.x.array.copy()


def figura_dupla(nome, desenhar, x_gap, titulo=None, cor_contorno="k"):
    """Figura com 2 painéis: domínio inteiro e zoom no dispositivo."""
    fig, axs = plt.subplots(1, 2, figsize=(11, 8),
                            gridspec_kw={"width_ratios": [1, 1.3]})
    rmax, zmin, zmax = viz.device_window(x_gap)
    for ax, janela in zip(axs, [(None, None, None), (rmax, zmin, zmax)]):
        desenhar(ax)
        viz.draw_device(ax, x_gap, color=cor_contorno)
        viz.set_window(ax, *janela)
    axs[0].set_title("Domínio completo")
    axs[1].set_title("Zoom no dispositivo")
    if titulo:
        fig.suptitle(titulo)
    viz.save(fig, OUT / nome)


def main():
    ap = argparse.ArgumentParser(description="Parte 1 do TC5")
    ap.add_argument("--x-gap-mm", type=float, default=cfg.X_GAP * 1e3)
    ap.add_argument("--current", type=float, default=1.0)
    ap.add_argument("--n-theta", type=int, default=361)
    args = ap.parse_args()
    x_gap = args.x_gap_mm * 1e-3
    I = args.current
    OUT.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------ simulação
    mesh, cell_tags, facet_tags = geometry.make_mesh(x_gap)
    n_cel = mesh.topology.index_map(2).size_local
    sol = solver.solve_magnetostatic(mesh, cell_tags, facet_tags, current=I)
    L = solver.inductance(sol)
    R = solver.ohmic_power(sol) / I ** 2

    # -------------------------------------------------------------- figuras
    pontos, tri, A_nod = dados_triangulacao(sol)
    t = viz.make_triangulation(pontos, tri)
    Bmag = solver.b_magnitude(sol)
    Hmag = np.linalg.norm(sol.H.x.array.reshape(-1, 2), axis=1)

    def desenha_malha(ax):
        viz.plot_mesh(ax, t, lw=0.12)
    figura_dupla("fig_malha.png", desenha_malha, x_gap,
                 titulo=f"Malha ({n_cel} elementos), x = {x_gap * 1e3:.1f} mm")

    def desenha_A(ax):
        viz.plot_nodal(ax, t, A_nod, "A' = r Aθ (Wb)")
        viz.plot_flux_lines(ax, t, A_nod, levels=25)
    figura_dupla("fig_A_fluxo.png", desenha_A, x_gap,
                 titulo=f"Potencial A' e linhas de fluxo (I = {I} A)")

    def desenha_B(ax):
        viz.plot_cells(ax, t, Bmag, "|B| (T)")
        viz.plot_flux_lines(ax, t, A_nod, levels=25, lw=0.4)
    figura_dupla("fig_B.png", desenha_B, x_gap,
                 titulo=f"Indução magnética |B| (I = {I} A)")

    def desenha_H(ax):
        viz.plot_cells(ax, t, Hmag, "|H| (A/m)", log=True)
    figura_dupla("fig_H.png", desenha_H, x_gap, cor_contorno="w",
                 titulo=f"Campo magnético |H| (I = {I} A)")

    # --------------------------------------------------------------- forças
    Fz_man, theta, Hr, Hz = postprocess.force_maxwell_manual(
        sol, x_gap, n_theta=args.n_theta)
    Fz_fen = postprocess.force_maxwell_fenics(sol)
    dif = 100 * (Fz_man - Fz_fen) / Fz_fen

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ang = np.degrees(theta)
    ax.plot(ang, Hr, label="Hr")
    ax.plot(ang, Hz, label="Hz")
    ax.set_xlabel("Ângulo θ na superfície da esfera (graus)\n"
                  "0° = equador, +90° = topo (lado do núcleo), -90° = base")
    ax.set_ylabel("Campo magnético (A/m)")
    ax.set_title(f"Hr e Hz na superfície externa da esfera (I = {I} A)")
    ax.grid(True)
    ax.legend()
    viz.save(fig, OUT / "fig_Hr_Hz_esfera.png")
    np.savetxt(OUT / "hr_hz_esfera.csv", np.column_stack([ang, Hr, Hz]),
               delimiter=",", header="theta_graus,Hr_A_m,Hz_A_m", comments="")

    # --------------------------------------- corrente de equilíbrio e o fio
    i_eq = postprocess.equilibrium_current(Fz_fen, I)
    i_eq_man = postprocess.equilibrium_current(Fz_man, I)
    fio = postprocess.wire_check(i_eq, R)

    # Refaz a simulação na corrente de equilíbrio e confere
    sol2 = solver.solve_magnetostatic(mesh, cell_tags, facet_tags, current=i_eq)
    Fz2_fen = postprocess.force_maxwell_fenics(sol2)
    Fz2_man = postprocess.force_maxwell_manual(sol2, x_gap, n_theta=args.n_theta)[0]
    b_nuc = postprocess.b_stats(sol2, cfg.TAG_CORE)
    b_esf = postprocess.b_stats(sol2, cfg.TAG_SPH)

    # ------------------------------------------------------------- resultados
    res = {
        "x_gap_mm": x_gap * 1e3,
        "n_celulas": int(n_cel),
        "corrente_referencia_A": I,
        "Fz_manual_N": Fz_man,
        "Fz_fenics_N": Fz_fen,
        "dif_manual_vs_fenics_pct": dif,
        "peso_N": cfg.PESO,
        "I_eq_A": i_eq,                    # usada na Parte 2
        "I_eq_manual_A": i_eq_man,
        "Fz_em_I_eq_fenics_N": Fz2_fen,
        "Fz_em_I_eq_manual_N": Fz2_man,
        "L_H": L,
        "R_ohm": R,
        "J_fio_A_mm2": fio["J_fio_A_mm2"],
        "P_W": fio["P_W"],
        "B_nucleo_max_T_em_I_eq": b_nuc["max"],
        "B_nucleo_med_T_em_I_eq": b_nuc["mean"],
        "B_esfera_max_T_em_I_eq": b_esf["max"],
        "B_esfera_med_T_em_I_eq": b_esf["mean"],
    }
    with open(OUT / "resultados.json", "w") as f:
        json.dump({k: float(v) for k, v in res.items()}, f, indent=2)

    linhas = [
        f"PARTE 1: x = {x_gap * 1e3:.2f} mm, malha com {n_cel} elementos",
        f"Fz com I = {I} A:",
        f"  Maxwell manual              = {Fz_man:.5e} N",
        f"  FacetNormal + assemble      = {Fz_fen:.5e} N",
        f"  diferença                   = {dif:+.2f} %",
        f"Peso da esfera                = {cfg.PESO:.5f} N",
        f"Corrente de equilíbrio        = {i_eq:.3f} A (assemble) | {i_eq_man:.3f} A (manual)",
        f"Conferência em I_eq           = {Fz2_fen:.5f} N (assemble) | {Fz2_man:.5f} N (manual)",
        f"L                             = {L * 1e3:.4f} mH",
        f"R da bobina                   = {R:.4f} ohm",
        f"Potência dissipada em I_eq    = {fio['P_W']:.2f} W",
        f"J no fio AWG {cfg.AWG_TIPO} em I_eq     = {fio['J_fio_A_mm2']:.2f} A/mm2",
        f"B em I_eq, núcleo: max = {b_nuc['max']:.4f} T | média = {b_nuc['mean']:.4f} T",
        f"B em I_eq, esfera: max = {b_esf['max']:.4f} T | média = {b_esf['mean']:.4f} T",
    ]
    texto = "\n".join(linhas)
    (OUT / "resumo.txt").write_text(texto + "\n", encoding="utf-8")
    print()
    print(texto)
    print(f"\nArquivos em: {OUT}")


if __name__ == "__main__":
    main()