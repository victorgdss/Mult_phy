"""
postprocess.py
Pós-processamento da solução magnetostática.

Força magnética na esfera (componente z) pelo tensor de Maxwell, em ar (mu0):

    Fz = ∮ mu0 [ (H.n) H_z - (1/2) |H|^2 n_z ] dS ,   dS = 2 pi r dl

onde n é a normal EXTERNA à esfera. Com o eixo z apontando para cima e o
núcleo acima da esfera, Fz > 0 significa atração (esfera sendo puxada).

Dois cálculos independentes:
    force_maxwell_manual : amostra Hr e Hz em pontos da superfície da esfera,
                           monta o integrando e integra (Simpson) em theta.
    force_maxwell_fenics : usa ufl.FacetNormal() e assemble_scalar() sobre as
                           facets do arco da esfera, pelo lado do ar (mesma
                           ideia da célula 55 do notebook da Aula 18).

Convenção do ângulo theta na superfície da esfera:
    theta = 0     equador da esfera (ponto de maior raio, r = R_SPH)
    theta = +90°  topo da esfera, ponto mais próximo do núcleo (sobre o eixo)
    theta = -90°  base da esfera (sobre o eixo)
    Ponto da superfície: r = R cos(theta), z = zc + R sin(theta), com zc o
    centro da esfera. Normal externa: n = (cos theta, sin theta).

Teste rápido (precisa do dolfinx):
    python postprocess.py                    # x = 5 mm, I = 1 A
    python postprocess.py --x-gap-mm 0.5
"""
import argparse

import numpy as np
import ufl
from scipy.integrate import simpson
from dolfinx import fem, geometry as dgeo

import config as cfg
import geometry
import solver


# ----------------------------------------------------------------------------
# Força: cálculo manual (tensor de Maxwell com Hr e Hz amostrados)
# ----------------------------------------------------------------------------
def _tag_por_celula(sol):
    """Vetor com a tag (1..4) de cada célula, indexado pelo índice da célula."""
    n_cel = sol.mesh.topology.index_map(sol.mesh.topology.dim).size_local
    tags = np.zeros(n_cel, dtype=np.int32)
    tags[sol.cell_tags.indices] = sol.cell_tags.values
    return tags


def sample_H_on_sphere(sol, x_gap, n_theta=181, delta=None):
    """
    Amostra Hr e Hz em pontos ligeiramente FORA da esfera (no ar).

    delta : afastamento radial a partir da superfície (m). Padrão: 0,5% de R_SPH.
            Deve ser menor que o tamanho dos elementos junto da esfera, para
            cair na primeira camada de células de ar.
    Devolve theta (rad), Hr, Hz (A/m).
    """
    if delta is None:
        delta = 0.005 * cfg.R_SPH
    zc = -x_gap - cfg.R_SPH
    theta = np.linspace(-np.pi / 2, np.pi / 2, n_theta)

    pts = np.zeros((n_theta, 3))
    pts[:, 0] = (cfg.R_SPH + delta) * np.cos(theta)
    pts[:, 1] = zc + (cfg.R_SPH + delta) * np.sin(theta)

    # Descobre em que célula cai cada ponto
    mesh = sol.mesh
    tree = dgeo.bb_tree(mesh, mesh.topology.dim)
    candidatos = dgeo.compute_collisions_points(tree, pts)
    colidindo = dgeo.compute_colliding_cells(mesh, candidatos, pts)
    tag_cel = _tag_por_celula(sol)

    celulas = []
    for i in range(n_theta):
        ar = [c for c in colidindo.links(i) if tag_cel[c] == cfg.TAG_AIR]
        if not ar:
            raise RuntimeError(
                f"Ponto {i} (theta = {np.degrees(theta[i]):.1f} graus) não caiu "
                f"em célula de ar. Reduza delta.")
        celulas.append(ar[0])

    H = sol.H.eval(pts, np.array(celulas, dtype=np.int32))
    return theta, H[:, 0], H[:, 1]


def force_maxwell_manual(sol, x_gap, n_theta=181, delta=None):
    """
    Fz (N) pelo tensor de Maxwell com integral manual em theta.
    Devolve (Fz, theta, Hr, Hz).
    """
    theta, Hr, Hz = sample_H_on_sphere(sol, x_gap, n_theta, delta)
    n_r, n_z = np.cos(theta), np.sin(theta)           # normal externa
    H2 = Hr ** 2 + Hz ** 2
    Hn = Hr * n_r + Hz * n_z
    r_s = cfg.R_SPH * np.cos(theta)                   # raio do ponto da superfície
    # dS = 2 pi r_s * (R dtheta)
    integrando = cfg.MU0 * (Hn * Hz - 0.5 * H2 * n_z) * 2 * np.pi * r_s * cfg.R_SPH
    Fz = simpson(integrando, x=theta)
    return Fz, theta, Hr, Hz


# ----------------------------------------------------------------------------
# Força: FacetNormal() + assemble_scalar()
# ----------------------------------------------------------------------------
def _medida_superficie_esfera(sol):
    """
    Medida ds que integra só nas facets do arco da esfera, usando em cada
    facet a célula de AR (assim a normal ufl.FacetNormal aponta do ar para a
    esfera; a normal externa à esfera é -n).
    """
    mesh = sol.mesh
    tdim = mesh.topology.dim
    fdim = tdim - 1
    mesh.topology.create_entities(fdim)
    mesh.topology.create_connectivity(fdim, tdim)
    mesh.topology.create_connectivity(tdim, fdim)
    f_para_c = mesh.topology.connectivity(fdim, tdim)
    c_para_f = mesh.topology.connectivity(tdim, fdim)

    tag_cel = _tag_por_celula(sol)
    facets = sol.facet_tags.find(cfg.TAG_SPH_SURF)
    if len(facets) == 0:
        raise RuntimeError("Nenhuma facet com a tag do arco da esfera.")

    entidades = []
    for f in facets:
        celulas = f_para_c.links(f)
        do_ar = [c for c in celulas if tag_cel[c] == cfg.TAG_AIR]
        assert len(do_ar) == 1, "facet do arco deveria ter exatamente 1 célula de ar"
        c = do_ar[0]
        idx_local = np.flatnonzero(c_para_f.links(c) == f)
        assert len(idx_local) == 1
        entidades.extend([c, idx_local[0]])

    entidades = np.array(entidades, dtype=np.int32)
    return ufl.Measure("ds", domain=mesh, subdomain_data=[(1, entidades)])


def force_maxwell_fenics(sol):
    """Fz (N) pelo tensor de Maxwell com FacetNormal() e assemble_scalar()."""
    ds2 = _medida_superficie_esfera(sol)
    n = ufl.FacetNormal(sol.mesh)          # aponta do ar para dentro da esfera
    n_ext = -n                             # normal externa à esfera
    H = sol.H
    integrando = (2 * np.pi * cfg.MU0 * sol.r
                  * (ufl.dot(H, n_ext) * H[1] - 0.5 * ufl.dot(H, H) * n_ext[1])
                  * ds2(1))
    return fem.assemble_scalar(fem.form(integrando))


# ----------------------------------------------------------------------------
# Indução magnética por região
# ----------------------------------------------------------------------------
def b_stats(sol, tag):
    """
    Bmax e Bmed (T) em uma região (cfg.TAG_CORE ou cfg.TAG_SPH).
    Bmax: maior |B| entre os elementos da região.
    Bmed: média PONDERADA PELO VOLUME, ∫ r |B| dr dz / ∫ r dr dz
          (a média simples por elemento dá peso demais às regiões com
          elementos pequenos).
    """
    Bmag = solver.b_magnitude(sol)
    bmax = Bmag[sol.cell_tags.find(tag)].max()
    num = fem.assemble_scalar(fem.form(
        sol.r * ufl.sqrt(ufl.dot(sol.B, sol.B)) * sol.dx(tag)))
    den = fem.assemble_scalar(fem.form(sol.r * sol.dx(tag)))
    return {"max": float(bmax), "mean": float(num / den)}


# ----------------------------------------------------------------------------
# Corrente de equilíbrio e verificação do fio
# ----------------------------------------------------------------------------
def equilibrium_current(fz, current):
    """
    Corrente em que Fz = peso. No problema LINEAR a força é proporcional a I^2,
    então I_eq = I * sqrt(P / Fz), com Fz calculado na corrente `current`.
    """
    return current * np.sqrt(cfg.PESO / abs(fz))


def wire_check(i_eq, resistance):
    """Densidade de corrente no fio e potência dissipada na corrente i_eq."""
    return {
        "J_fio_A_mm2": i_eq / (cfg.fio.area * 1e6),
        "P_W": resistance * i_eq ** 2,
    }


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Teste rápido do pós-processamento")
    ap.add_argument("--x-gap-mm", type=float, default=cfg.X_GAP * 1e3)
    ap.add_argument("--current", type=float, default=1.0)
    ap.add_argument("--n-theta", type=int, default=181)
    args = ap.parse_args()
    x_gap = args.x_gap_mm * 1e-3
    I = args.current

    mesh, cell_tags, facet_tags = geometry.make_mesh(x_gap)
    sol = solver.solve_magnetostatic(mesh, cell_tags, facet_tags, current=I)

    Fz_man, theta, Hr, Hz = force_maxwell_manual(sol, x_gap, n_theta=args.n_theta)
    Fz_fen = force_maxwell_fenics(sol)
    dif = 100 * (Fz_man - Fz_fen) / Fz_fen

    print(f"x = {x_gap * 1e3:.3f} mm | I = {I} A | {mesh.topology.index_map(2).size_local} células")
    print(f"  Fz (Maxwell manual)         = {Fz_man:.5e} N")
    print(f"  Fz (FacetNormal + assemble) = {Fz_fen:.5e} N")
    print(f"  diferença manual vs assemble = {dif:+.2f} %")
    print(f"  Hz no topo da esfera (theta = +90°) = {Hz[-1]:.1f} A/m")

    for nome, tag in (("núcleo", cfg.TAG_CORE), ("esfera", cfg.TAG_SPH)):
        s = b_stats(sol, tag)
        print(f"  B {nome}: max = {s['max']:.4f} T | média (ponderada) = {s['mean']:.4f} T")

    # Corrente de equilíbrio (a partir de cada valor de força) e checagem
    for nome, F in (("manual", Fz_man), ("assemble", Fz_fen)):
        print(f"  I_eq a partir do Fz {nome:<8} = {equilibrium_current(F, I):.3f} A")

    i_eq = equilibrium_current(Fz_fen, I)
    R = solver.ohmic_power(sol) / I ** 2
    w = wire_check(i_eq, R)
    print(f"  Fio AWG {cfg.AWG_TIPO} em I_eq = {i_eq:.2f} A: "
          f"J = {w['J_fio_A_mm2']:.2f} A/mm2 | P = {w['P_W']:.2f} W")

    # Verificação: refaz a simulação em I_eq e confere se a força = peso
    sol2 = solver.solve_magnetostatic(mesh, cell_tags, facet_tags, current=i_eq)
    F2 = force_maxwell_fenics(sol2)
    print(f"  Conferência: Fz(I_eq) = {F2:.5f} N | peso = {cfg.PESO:.5f} N")