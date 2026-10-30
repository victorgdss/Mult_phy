"""
parte2_parametros.py
Parte 2 do TC5: parâmetros concentrados do levitador em função da posição x.

Para cada posição x (10 pontos linearmente espaçados de 0,1 x_eq a 3,0 x_eq,
mais o próprio ponto de equilíbrio x_eq), com a corrente de equilíbrio da
Parte 1 (modelo linear, então L e Fm não dependem de como a corrente é escolhida
além do fator I^2 na força):

    L(x)    indutância, pela energia magnética
    Fm(x)   força magnética (FacetNormal + assemble); a do tensor de Maxwell
            manual também é calculada, só para checar a diferença
    Bl(x)   fator de força, Bl = Fm / I   (definição da Aula 18)
    Bmax, Bmed no núcleo e na esfera

Comparação com o modelo planar da Aula 18:
    L(x)  = k / x
    Fm    = (k/2) (I/x)^2
    k_L   = x L(x)                 k_F = 2 Fm (x/I)^2
    k_th  = mu0 N^2 S / 2          (forma fechada, com S = área da seção do núcleo)

Saídas em results/parte2/: tabelas (CSV e texto), figuras (PNG).

Uso:
    python parte2_parametros.py                 # lê I_eq de results/parte1
    python parte2_parametros.py --current 6.23  # corrente informada na mão
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import config as cfg
import geometry
import solver
import postprocess
import viz

BASE = Path(__file__).resolve().parent / "results"
OUT = BASE / "parte2"


def carregar_corrente(valor_cli):
    if valor_cli is not None:
        return valor_cli
    caminho = BASE / "parte1" / "resultados.json"
    if not caminho.exists():
        sys.exit("Não achei results/parte1/resultados.json. Rode parte1_base.py "
                 "primeiro ou passe --current.")
    with open(caminho) as f:
        return json.load(f)["I_eq_A"]


def posicoes(inicio, fim, n, incluir_equilibrio=True):
    """Posições x (m): n pontos entre inicio*X_GAP e fim*X_GAP, mais x_eq."""
    xs = list(np.linspace(inicio * cfg.X_GAP, fim * cfg.X_GAP, n))
    if incluir_equilibrio and not any(abs(x - cfg.X_GAP) < 1e-9 for x in xs):
        xs.append(cfg.X_GAP)
    return np.array(sorted(xs))


def calcular_ponto(x, I, n_theta):
    """Resolve o problema em x e devolve um dicionário com os parâmetros."""
    mesh, cell_tags, facet_tags = geometry.make_mesh(x)
    sol = solver.solve_magnetostatic(mesh, cell_tags, facet_tags, current=I)
    Fm = postprocess.force_maxwell_fenics(sol)
    Fm_man = postprocess.force_maxwell_manual(sol, x, n_theta=n_theta)[0]
    b_nuc = postprocess.b_stats(sol, cfg.TAG_CORE)
    b_esf = postprocess.b_stats(sol, cfg.TAG_SPH)
    return {
        "x_m": x,
        "n_celulas": mesh.topology.index_map(2).size_local,
        "L_H": solver.inductance(sol),
        "Fm_N": Fm,
        "Fm_manual_N": Fm_man,
        "Bl_N_A": Fm / I,
        "Bmax_nucleo_T": b_nuc["max"], "Bmed_nucleo_T": b_nuc["mean"],
        "Bmax_esfera_T": b_esf["max"], "Bmed_esfera_T": b_esf["mean"],
    }


def pos_processar(linhas, I):
    """Acrescenta k_L, k_F, razão, expoente local da força e diferença manual."""
    x = np.array([l["x_m"] for l in linhas])
    L = np.array([l["L_H"] for l in linhas])
    Fm = np.array([l["Fm_N"] for l in linhas])
    kL = x * L
    kF = 2 * Fm * (x / I) ** 2
    # expoente local n em Fm ~ x^(-n): o modelo planar prevê n = 2
    n_loc = -np.gradient(np.log(Fm), np.log(x))
    for i, l in enumerate(linhas):
        l["k_L"], l["k_F"], l["razao_kL_kF"] = kL[i], kF[i], kL[i] / kF[i]
        l["expoente_F"] = n_loc[i]
        l["dif_manual_pct"] = 100 * (l["Fm_manual_N"] - l["Fm_N"]) / l["Fm_N"]

    # Ajuste L(x) = L0 + a/x (indutância de dispersão L0 + parcela que acopla na esfera)
    M = np.column_stack([np.ones_like(x), 1 / x])
    (L0, a), *_ = np.linalg.lstsq(M, L, rcond=None)
    erro_ajuste = np.max(np.abs(M @ np.array([L0, a]) - L) / L)
    return {"L0_H": L0, "a_H_m": a, "erro_max_ajuste": erro_ajuste}


def salvar_csv(linhas, caminho):
    colunas = ["x_m", "L_H", "Fm_N", "Fm_manual_N", "dif_manual_pct", "Bl_N_A",
               "Bmax_nucleo_T", "Bmed_nucleo_T", "Bmax_esfera_T", "Bmed_esfera_T",
               "k_L", "k_F", "razao_kL_kF", "expoente_F", "n_celulas"]
    with open(caminho, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(colunas)
        for l in linhas:
            w.writerow([f"{l[c]:.8e}" if c != "n_celulas" else int(l[c]) for c in colunas])


def tabelas_texto(linhas, eq):
    t = []
    t.append("TABELA 1: parâmetros concentrados (* = posição de equilíbrio)")
    t.append(f"{'x (mm)':>9}{'L (mH)':>10}{'Fm (N)':>11}{'Bl (N/A)':>11}{'dif manual %':>14}")
    for l, e in zip(linhas, eq):
        t.append(f"{l['x_m'] * 1e3:>8.3f}{'*' if e else ' '}{l['L_H'] * 1e3:>10.4f}"
                 f"{l['Fm_N']:>11.5f}{l['Bl_N_A']:>11.5f}{l['dif_manual_pct']:>14.2f}")
    t.append("")
    t.append("TABELA 2: indução magnética (T)")
    t.append(f"{'x (mm)':>9}{'Bmax núcleo':>14}{'Bmed núcleo':>14}{'Bmax esfera':>14}{'Bmed esfera':>14}")
    for l, e in zip(linhas, eq):
        t.append(f"{l['x_m'] * 1e3:>8.3f}{'*' if e else ' '}{l['Bmax_nucleo_T']:>14.4f}"
                 f"{l['Bmed_nucleo_T']:>14.4f}{l['Bmax_esfera_T']:>14.4f}{l['Bmed_esfera_T']:>14.4f}")
    t.append("")
    t.append("TABELA 3: constante k do modelo planar, k = xL e k = 2 Fm (x/I)^2")
    t.append(f"{'x (mm)':>9}{'k_L':>13}{'k_F':>13}{'k_L/k_F':>10}{'expoente F':>12}")
    for l, e in zip(linhas, eq):
        t.append(f"{l['x_m'] * 1e3:>8.3f}{'*' if e else ' '}{l['k_L']:>13.3e}{l['k_F']:>13.3e}"
                 f"{l['razao_kL_kF']:>10.1f}{l['expoente_F']:>12.2f}")
    return "\n".join(t)


def figuras(linhas, I, k_th):
    x = np.array([l["x_m"] for l in linhas])
    xmm = x * 1e3
    L = np.array([l["L_H"] for l in linhas])
    Fm = np.array([l["Fm_N"] for l in linhas])
    kL = np.array([l["k_L"] for l in linhas])
    kF = np.array([l["k_F"] for l in linhas])
    eq = np.isclose(x, cfg.X_GAP, rtol=0, atol=1e-9)
    xf = np.linspace(x.min(), x.max(), 300)

    # --- L(x)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(xmm, L * 1e3, "o-", label="FEM")
    ax.loglog(xf * 1e3, k_th / xf * 1e3, "--", label="planar analítico, L = k/x")
    ax.loglog(xmm[eq], L[eq] * 1e3, "r*", ms=12, label="equilíbrio")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("L (mH)")
    ax.set_title("Indutância em função da posição"); ax.grid(True, which="both"); ax.legend()
    viz.save(fig, OUT / "fig_L_x.png")

    # --- Fm(x)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(xmm, Fm, "o-", label="FEM")
    ax.loglog(xf * 1e3, k_th / 2 * (I / xf) ** 2, "--", label="planar analítico, Fm = (k/2)(I/x)²")
    ax.axhline(cfg.PESO, color="gray", ls=":", label="peso da esfera")
    ax.loglog(xmm[eq], Fm[eq], "r*", ms=12, label="equilíbrio")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("Fm (N)")
    ax.set_title(f"Força magnética em função da posição (I = {I:.2f} A)")
    ax.grid(True, which="both"); ax.legend()
    viz.save(fig, OUT / "fig_Fm_x.png")

    # --- k
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.semilogy(xmm, kL, "o-", label="k = x L(x)")
    ax.semilogy(xmm, kF, "s-", label="k = 2 Fm (x/I)²")
    ax.axhline(k_th, color="gray", ls="--", label="k = μ0 N² S / 2")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("k (H·m)")
    ax.set_title("Constante k equivalente do modelo planar")
    ax.grid(True, which="both"); ax.legend()
    viz.save(fig, OUT / "fig_k_x.png")

    # --- B
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, regiao, rot in zip(axs, ("nucleo", "esfera"), ("núcleo", "esfera")):
        ax.plot(xmm, [l[f"Bmax_{regiao}_T"] for l in linhas], "o-", label="Bmax")
        ax.plot(xmm, [l[f"Bmed_{regiao}_T"] for l in linhas], "s-", label="Bmed")
        ax.set_xlabel("x (mm)"); ax.set_ylabel("B (T)")
        ax.set_title(f"Indução no {rot}"); ax.grid(True); ax.legend()
    viz.save(fig, OUT / "fig_B_x.png")


def main():
    ap = argparse.ArgumentParser(description="Parte 2 do TC5")
    ap.add_argument("--current", type=float, default=None,
                    help="corrente (A); padrão: I_eq de results/parte1")
    ap.add_argument("--inicio", type=float, default=0.1, help="início, em múltiplos de x")
    ap.add_argument("--fim", type=float, default=3.0, help="fim, em múltiplos de x")
    ap.add_argument("--n-pontos", type=int, default=10)
    ap.add_argument("--n-theta", type=int, default=361)
    args = ap.parse_args()

    I = carregar_corrente(args.current)
    OUT.mkdir(parents=True, exist_ok=True)
    xs = posicoes(args.inicio, args.fim, args.n_pontos)
    print(f"Corrente usada: I = {I:.4f} A | {len(xs)} posições "
          f"({xs.min() * 1e3:.2f} a {xs.max() * 1e3:.2f} mm)")

    linhas = []
    for i, x in enumerate(xs, 1):
        l = calcular_ponto(x, I, args.n_theta)
        linhas.append(l)
        print(f"  [{i:2d}/{len(xs)}] x = {x * 1e3:6.3f} mm | {l['n_celulas']:6d} células | "
              f"L = {l['L_H'] * 1e3:.4f} mH | Fm = {l['Fm_N']:.5f} N")

    ajuste = pos_processar(linhas, I)
    eq = [abs(l["x_m"] - cfg.X_GAP) < 1e-9 for l in linhas]
    S = np.pi * cfg.R_CORE ** 2
    k_th = cfg.MU0 * cfg.N_TURNS ** 2 * S / 2

    salvar_csv(linhas, OUT / "dados_parte2.csv")
    figuras(linhas, I, k_th)

    texto = tabelas_texto(linhas, eq)
    i_eq_row = linhas[int(np.argmax(eq))] if any(eq) else None
    extra = [
        "",
        f"k teórico do planar (mu0 N^2 S / 2, S = pi r_core^2) = {k_th:.3e} H.m",
        f"Ajuste L(x) = L0 + a/x : L0 = {ajuste['L0_H'] * 1e3:.4f} mH, "
        f"a = {ajuste['a_H_m']:.3e} H.m (erro máx. do ajuste {100 * ajuste['erro_max_ajuste']:.2f} %)",
    ]
    if i_eq_row:
        extra.append(f"Conferência no equilíbrio: Fm = {i_eq_row['Fm_N']:.5f} N | "
                     f"peso = {cfg.PESO:.5f} N")
    texto += "\n" + "\n".join(extra)
    (OUT / "resumo.txt").write_text(texto + "\n", encoding="utf-8")
    print()
    print(texto)
    print(f"\nArquivos em: {OUT}")


if __name__ == "__main__":
    main()
