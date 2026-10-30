"""
viz.py
Funções de gráfico em matplotlib puro (não depende do FEniCS).

Convenção: os dados chegam em metros e os gráficos saem em milímetros.
Eixo horizontal = r, eixo vertical = z. Só o lado r >= 0 é simulado.
O backend é o "Agg": os scripts salvam PNG e não abrem janela.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Rectangle, Arc

import config as cfg

MM = 1e3


def make_triangulation(points, triangles):
    """points: (N,2) em metros; triangles: (M,3) índices dos vértices."""
    return mtri.Triangulation(points[:, 0] * MM, points[:, 1] * MM, triangles)


def draw_device(ax, x_gap, color="k", lw=0.9):
    """Contorno do núcleo, da bobina (tracejado) e da esfera."""
    ax.add_patch(Rectangle((0, 0), cfg.R_CORE * MM, cfg.Z_CORE * MM,
                           fill=False, ec=color, lw=lw))
    ax.add_patch(Rectangle((cfg.R_CORE * MM, cfg.COIL_Z0 * MM),
                           cfg.COIL_THICKNESS * MM, cfg.Z_COIL * MM,
                           fill=False, ec=color, lw=lw, ls="--"))
    zc = (-x_gap - cfg.R_SPH) * MM
    d = 2 * cfg.R_SPH * MM
    ax.add_patch(Arc((0, zc), d, d, theta1=-90, theta2=90, ec=color, lw=lw))


def set_window(ax, rmax_mm=None, zmin_mm=None, zmax_mm=None):
    """Recorta a janela do gráfico (None = deixa automático)."""
    ax.set_xlim(0, rmax_mm if rmax_mm is not None else cfg.R_AIR * MM)
    ax.set_ylim(zmin_mm if zmin_mm is not None else cfg.Z_AIR_MIN * MM,
                zmax_mm if zmax_mm is not None else cfg.Z_AIR_MAX * MM)
    ax.set_aspect("equal")
    ax.set_xlabel("r (mm)")
    ax.set_ylabel("z (mm)")


def device_window(x_gap, margem_mm=6.0):
    """Janela que enquadra o dispositivo (núcleo, bobina e esfera)."""
    rmax = max(cfg.R_COIL, cfg.R_SPH) * MM + margem_mm
    zmin = (-x_gap - 2 * cfg.R_SPH) * MM - margem_mm
    zmax = cfg.Z_CORE * MM + margem_mm
    return rmax, zmin, zmax


def plot_mesh(ax, tri, lw=0.15):
    ax.triplot(tri, lw=lw, color="k")


def plot_cells(ax, tri, values, label, log=False, cmap="viridis"):
    """Campo constante por elemento (DG0). log=True mostra log10 do valor."""
    v = np.asarray(values, dtype=float)
    if log:
        v = np.log10(np.clip(v, 1e-12, None))
        label = f"log10({label})"
    pc = ax.tripcolor(tri, facecolors=v, cmap=cmap, edgecolors="none")
    cb = ax.figure.colorbar(pc, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label(label)
    return pc


def plot_nodal(ax, tri, values, label, cmap="viridis"):
    """Campo nodal (P1), interpolado suavemente."""
    pc = ax.tripcolor(tri, np.asarray(values, dtype=float), cmap=cmap, shading="gouraud")
    cb = ax.figure.colorbar(pc, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label(label)
    return pc


def plot_flux_lines(ax, tri, a_nodal, levels=25, color="w", lw=0.6):
    """Linhas de nível de A' = linhas de fluxo de B."""
    ax.tricontour(tri, np.asarray(a_nodal, dtype=float), levels=levels,
                  colors=color, linewidths=lw)


def save(fig, caminho, dpi=200):
    fig.tight_layout()
    fig.savefig(caminho, dpi=dpi)
    plt.close(fig)
    print("figura salva:", caminho)