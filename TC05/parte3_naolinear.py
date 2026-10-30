"""
parte3_naolinear.py
Parte 3 do TC5: efeitos da não linearidade (análise sobre a curva B x H).

Entrada: results/parte2/dados_parte2.csv (Bmax e Bmed no núcleo e na esfera,
em cada posição, da simulação LINEAR com mur = 1200, na corrente de equilíbrio).
Este script não resolve FEM: ele só usa esses B e a curva da Tabela 3.

Para cada B (4 valores por posição) o script determina:
    H_linear  = B / (mur_linear mu0)   H que o FEM linear realmente usa
                                       (ponto sobre a reta B = 1200 mu0 H)
    H_curva   = H(B) da curva real     ponto de operação na curva do ferro puro
    mur_real  = B / (mu0 H_curva)      permeabilidade secante esperada no
                                       problema não linear
    razão     = mur_real / mur_linear

Saídas em results/parte3/: pontos_operacao.csv, resumo.txt, fig_BH_pontos.png,
fig_mur.png.

Uso:
    python parte3_naolinear.py
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import config as cfg
import bh_curve as bh
import viz

BASE = Path(__file__).resolve().parent / "results"
OUT = BASE / "parte3"

# (coluna do CSV da Parte 2, rótulo, região, cor, marcador)
CATEGORIAS = [
    ("Bmax_nucleo_T", "Bmax núcleo", "nucleo", "tab:blue", "o"),
    ("Bmed_nucleo_T", "Bmed núcleo", "nucleo", "tab:cyan", "s"),
    ("Bmax_esfera_T", "Bmax esfera", "esfera", "tab:red", "^"),
    ("Bmed_esfera_T", "Bmed esfera", "esfera", "tab:orange", "v"),
]
MUR_LINEAR = {"nucleo": cfg.MUR_CORE, "esfera": cfg.MUR_SPH}


def carregar_parte2(caminho):
    if not Path(caminho).exists():
        sys.exit(f"Não achei {caminho}. Rode parte2_parametros.py primeiro.")
    with open(caminho) as f:
        linhas = list(csv.DictReader(f))
    return {k: np.array([float(l[k]) for l in linhas]) for k in linhas[0]}


def corrente_equilibrio():
    caminho = BASE / "parte1" / "resultados.json"
    if caminho.exists():
        with open(caminho) as f:
            return json.load(f)["I_eq_A"]
    return None


def calcular(dados):
    """Dicionário: rótulo -> arrays B, H_linear, H_curva, mur_real, razão."""
    res = {}
    for col, rotulo, regiao, _, _ in CATEGORIAS:
        B = dados[col]
        mur_lin = MUR_LINEAR[regiao]
        H_lin = B / (mur_lin * cfg.MU0)
        H_cur = bh.h_from_b(B)
        mur_real = bh.mur_secante(B)
        res[rotulo] = dict(B=B, H_lin=H_lin, H_cur=H_cur, mur_real=mur_real,
                           razao=mur_real / mur_lin, mur_lin=mur_lin)
    return res


def salvar_csv(x, res, caminho):
    with open(caminho, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "grandeza", "B_T", "H_linear_A_m", "H_curva_A_m",
                    "mur_linear", "mur_real_secante", "razao"])
        for i, xi in enumerate(x):
            for rotulo, r in res.items():
                w.writerow([f"{xi * 1e3:.4f}", rotulo, f"{r['B'][i]:.6f}",
                            f"{r['H_lin'][i]:.4f}", f"{r['H_cur'][i]:.4f}",
                            r["mur_lin"], f"{r['mur_real'][i]:.1f}", f"{r['razao'][i]:.2f}"])


def tabela_texto(x, res):
    rot = list(res)
    cab = f"{'x (mm)':>9}" + "".join(f"{r:>13}{'mur':>8}" for r in rot)
    linhas = ["TABELA: B da simulação linear e mur esperado na curva real (secante)", cab]
    for i, xi in enumerate(x):
        marca = "*" if abs(xi - cfg.X_GAP) < 1e-9 else " "
        linhas.append(f"{xi * 1e3:>8.3f}{marca}" + "".join(
            f"{res[r]['B'][i]:>12.4f}T{res[r]['mur_real'][i]:>8.0f}" for r in rot))
    return "\n".join(linhas)


def figura_bh(x, res):
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))

    # Painel 1: curva inteira, pontos de operação na curva real
    ax = axs[0]
    ax.plot(bh.H_TAB, bh.B_TAB, "-", color="gray", label="curva do ferro puro (Tab. 3)")
    for _, rotulo, _, cor, mk in CATEGORIAS:
        ax.plot(res[rotulo]["H_cur"], res[rotulo]["B"], mk, color=cor, ms=5, label=rotulo)
    ax.set_xlabel("H (A/m)"); ax.set_ylabel("B (T)")
    ax.set_title("Curva B×H completa e pontos de operação")
    ax.grid(True); ax.legend(loc="lower right", fontsize=8)

    # Painel 2: zoom em baixo H, com a reta B = 1200 mu0 H usada no FEM linear
    ax = axs[1]
    Bf = np.linspace(0, 1.5, 600)
    ax.plot(bh.h_from_b(Bf), Bf, "-", color="gray", lw=2, label="curva real (Tab. 3)")
    Hl = np.linspace(0, 400, 50)
    ax.plot(Hl, cfg.MUR_CORE * cfg.MU0 * Hl, "k--", lw=1.2,
            label=f"FEM linear: B = {cfg.MUR_CORE}·μ0·H")
    for _, rotulo, _, cor, mk in CATEGORIAS:
        r = res[rotulo]
        ax.plot(r["H_lin"], r["B"], mk, mfc="none", color=cor, ms=6)
        ax.plot(r["H_cur"], r["B"], mk, color=cor, ms=5, label=f"{rotulo}")
    ax.plot([], [], "ko", mfc="none", label="ponto do FEM linear (vazado)")
    ax.plot([], [], "ko", label="mesmo B na curva real (cheio)")
    ax.set_xlim(0, 400); ax.set_ylim(0, 0.65)
    ax.set_xlabel("H (A/m)"); ax.set_ylabel("B (T)")
    ax.set_title("Zoom: mesmo B, H diferente")
    ax.grid(True); ax.legend(loc="lower right", fontsize=7)
    viz.save(fig, OUT / "fig_BH_pontos.png")


def figura_mur(x, res):
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axs[0]
    Bf = np.linspace(0.005, bh.B_TAB[-1], 800)
    ax.semilogx(Bf, bh.mur_secante(Bf), "-", color="gray", label="mur secante da curva real")
    ax.axhline(cfg.MUR_CORE, color="k", ls="--", label=f"mur linear = {cfg.MUR_CORE}")
    for _, rotulo, _, cor, mk in CATEGORIAS:
        ax.semilogx(res[rotulo]["B"], res[rotulo]["mur_real"], mk, color=cor, ms=5, label=rotulo)
    ax.axvline(bh.b_joelho(), color="r", ls=":", label="mur cai à metade")
    ax.set_xlabel("B (T)"); ax.set_ylabel("permeabilidade relativa")
    ax.set_title("Permeabilidade secante em função de B")
    ax.grid(True, which="both"); ax.legend(fontsize=8)

    ax = axs[1]
    xmm = x * 1e3
    for _, rotulo, _, cor, mk in CATEGORIAS:
        ax.plot(xmm, res[rotulo]["mur_real"], mk + "-", color=cor, ms=5, label=rotulo)
    ax.axhline(cfg.MUR_CORE, color="k", ls="--", label=f"mur linear = {cfg.MUR_CORE}")
    ax.set_ylim(0, 1.15 * max(r["mur_real"].max() for r in res.values()))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("mur esperado (curva real)")
    ax.set_title("Permeabilidade esperada nos pontos de operação")
    ax.grid(True); ax.legend(fontsize=8)
    viz.save(fig, OUT / "fig_mur.png")


def main():
    ap = argparse.ArgumentParser(description="Parte 3 do TC5")
    ap.add_argument("--parte2-csv", default=str(BASE / "parte2" / "dados_parte2.csv"))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    dados = carregar_parte2(args.parte2_csv)
    x = dados["x_m"]
    res = calcular(dados)
    salvar_csv(x, res, OUT / "pontos_operacao.csv")
    figura_bh(x, res)
    figura_mur(x, res)

    # ---------------------------------------------------------------- resumo
    todos_B = np.concatenate([r["B"] for r in res.values()])
    todos_mur = np.concatenate([r["mur_real"] for r in res.values()])
    todos_razao = np.concatenate([r["razao"] for r in res.values()])
    b_nuc_max = dados["Bmax_nucleo_T"].max()
    b_j = bh.b_joelho()
    i_eq = corrente_equilibrio()

    linhas = [tabela_texto(x, res), ""]
    linhas += [
        f"B nos pontos de operação: de {todos_B.min():.4f} T a {todos_B.max():.4f} T",
        f"mur esperado (curva real): de {todos_mur.min():.0f} a {todos_mur.max():.0f} "
        f"(mur inicial da curva = {bh.MUR_INICIAL:.0f})",
        f"mur linear usado no FEM = {cfg.MUR_CORE}; razão real/linear: "
        f"de {todos_razao.min():.1f} a {todos_razao.max():.1f}",
        f"H que o FEM linear usa vs H da curva real, no pior caso (Bmax núcleo = {b_nuc_max:.3f} T): "
        f"{b_nuc_max / (cfg.MUR_CORE * cfg.MU0):.0f} A/m vs {float(bh.h_from_b(b_nuc_max)):.1f} A/m",
        f"B em que mur secante cai à metade do inicial (joelho de saturação): {b_j:.2f} T",
        f"Maior B de toda a simulação / joelho = {todos_B.max() / b_j:.2f}",
    ]
    if i_eq:
        linhas.append(
            f"Estimativa (modelo linear, B proporcional a I): o núcleo só chegaria ao joelho com "
            f"I ~ {i_eq * b_j / b_nuc_max:.1f} A (I_eq = {i_eq:.2f} A)")
    texto = "\n".join(linhas)
    (OUT / "resumo.txt").write_text(texto + "\n", encoding="utf-8")
    print(texto)
    print(f"\nArquivos em: {OUT}")


if __name__ == "__main__":
    main()
