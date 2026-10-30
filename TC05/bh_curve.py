"""
bh_curve.py
Curva de magnetização do ferro puro (Tabela 3 do enunciado, biblioteca do FEMM)
e funções de interpolação. Não depende do FEniCS.

Definições usadas:
    permeabilidade relativa SECANTE   mur(B) = B / (mu0 H(B))
    (é a que entra em B = mu0 mur H e a que um solver não linear usa a cada
    iteração; a permeabilidade diferencial dB/dH é outra coisa)
"""
import numpy as np

import config as cfg

# Tabela 3 do enunciado
B_TAB = np.array([
    0.000000, 0.227065, 0.454130, 0.681195, 0.908260, 1.135330, 1.362390,
    1.589350, 1.812360, 2.010040, 2.133160, 2.199990, 2.254790, 2.299930,
    2.342510, 2.378760, 2.415010, 2.451260, 2.487500, 2.523750, 2.560000])   # (T)
H_TAB = np.array([
    0.0000, 13.8984, 27.7967, 42.3974, 61.4157, 82.3824, 144.6690,
    897.7600, 4581.7400, 17736.0000, 41339.3000, 68321.8000, 95685.5000,
    123355.0000, 151083.0000, 178954.0000, 206825.0000, 234696.0000,
    262568.0000, 290439.0000, 318310.0000])                                   # (A/m)

# Permeabilidade relativa no início da curva (primeiro segmento)
MUR_INICIAL = B_TAB[1] / (cfg.MU0 * H_TAB[1])


def h_from_b(B):
    """H (A/m) correspondente a B (T) na curva. Acima da tabela, extrapola com
    a inclinação do vácuo (material totalmente saturado)."""
    B = np.abs(np.asarray(B, dtype=float))
    H = np.interp(B, B_TAB, H_TAB)
    return np.where(B > B_TAB[-1], H_TAB[-1] + (B - B_TAB[-1]) / cfg.MU0, H)


def mur_secante(B):
    """Permeabilidade relativa secante B / (mu0 H(B)). Em B ~ 0 devolve o valor inicial."""
    B = np.abs(np.asarray(B, dtype=float))
    H = h_from_b(B)
    with np.errstate(divide="ignore", invalid="ignore"):
        mur = B / (cfg.MU0 * H)
    return np.where(B > 1e-12, mur, MUR_INICIAL)


def b_joelho(fracao=0.5):
    """Menor B em que mur secante cai abaixo de `fracao` do valor inicial."""
    B = np.linspace(1e-3, B_TAB[-1], 20000)
    abaixo = mur_secante(B) < fracao * MUR_INICIAL
    return float(B[np.argmax(abaixo)])
