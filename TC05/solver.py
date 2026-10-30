"""
solver.py
Problema magnetostático axissimétrico do levitador, resolvido em A' = r * A_theta.

Equação no plano (r, z), com mu = mu0 * mur:

    -div( (1 / (mu r)) grad A' ) = J_theta            (forma forte)

Forma fraca (w = função teste, A' = 0 no contorno do domínio):

    ∫ (1 / (mu r)) grad(w) . grad(A') dr dz = ∫ w J dr dz

Campos obtidos a partir de A':

    B_r = -(1/r) dA'/dz        B_z = (1/r) dA'/dr        H = B / mu

Este arquivo tem:
    solve_magnetostatic(...)  resolve o problema para uma malha e uma corrente
    magnetic_energy / inductance / ohmic_power   grandezas globais

Teste rápido (gera a malha, resolve e imprime valores de conferência):
    python solver.py                       # x = 5 mm, I = 1 A
    python solver.py --x-gap-mm 0.5 --current 6
"""
import argparse
from dataclasses import dataclass
from typing import Any

import numpy as np
import ufl
from petsc4py import PETSc
from dolfinx import fem, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from dolfinx.mesh import locate_entities_boundary

import config as cfg
import geometry


@dataclass
class MagnetoSolution:
    """Tudo que os pós-processamentos precisam."""
    mesh: Any
    cell_tags: Any
    facet_tags: Any
    current: float
    r: Any          # coordenada radial (expressão UFL)
    dx: Any         # medida de volume com as tags das regiões
    A: Any          # A' = r A_theta   (Lagrange P2)
    mu: Any         # permeabilidade absoluta (DG0)
    J: Any          # densidade de corrente (DG0)
    B: Any          # indução magnética (DG0 vetorial: Br, Bz)
    H: Any          # campo magnético (DG0 vetorial: Hr, Hz)


def solve_magnetostatic(mesh, cell_tags, facet_tags=None, current=1.0,
                        mur_core=cfg.MUR_CORE, mur_sph=cfg.MUR_SPH):
    """Resolve o problema magnetostático linear para a corrente dada (A)."""
    tdim = mesh.topology.dim

    # Espaço de funções: Lagrange P2 para A'
    V = fem.functionspace(mesh, ("Lagrange", 2))
    x = ufl.SpatialCoordinate(mesh)
    r = x[0]

    # Condição de contorno de Dirichlet: A' = 0 em TODA a fronteira da caixa
    # (inclui o eixo r = 0, onde A' = r A_theta = 0 de qualquer forma)
    facets = locate_entities_boundary(mesh, tdim - 1,
                                      lambda x: np.full(x.shape[1], True))
    dofs = fem.locate_dofs_topological(V, tdim - 1, facets)
    bc = fem.dirichletbc(default_scalar_type(0), dofs, V)

    # Permeabilidade absoluta, constante por elemento (DG0).
    # Em DG0 o índice do grau de liberdade é o índice da célula (como no
    # notebook do professor).
    Q = fem.functionspace(mesh, ("DG", 0))
    mu = fem.Function(Q)
    mu.x.array[:] = cfg.MU0                                  # ar e bobina
    mu.x.array[cell_tags.find(cfg.TAG_CORE)] = mur_core * cfg.MU0
    mu.x.array[cell_tags.find(cfg.TAG_SPH)] = mur_sph * cfg.MU0

    dx = ufl.Measure("dx", domain=mesh, subdomain_data=cell_tags)

    # Densidade de corrente na bobina: J = N I / (área da bobina no plano r,z)
    area_coil = fem.assemble_scalar(fem.form(1 * dx(cfg.TAG_COIL)))
    J = fem.Function(Q)
    J.x.array[:] = 0.0
    J.x.array[cell_tags.find(cfg.TAG_COIL)] = cfg.N_TURNS * current / area_coil

    # Forma fraca
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = (1 / mu) * (1 / r) * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx
    L = J * v * dx

    # Solver direto (LU): o problema é pequeno e o contraste de mu é grande
    problem = LinearProblem(a, L, bcs=[bc],
                            petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    A = problem.solve()
    A.name = "A_prime"

    # B = rot(A) e H = B / mu, avaliados no centro de cada elemento (DG0)
    W = fem.functionspace(mesh, ("DG", 0, (mesh.geometry.dim,)))
    B = fem.Function(W)
    B.interpolate(fem.Expression(
        ufl.as_vector((-(1 / r) * A.dx(1), (1 / r) * A.dx(0))),
        W.element.interpolation_points()))
    H = fem.Function(W)
    H.interpolate(fem.Expression(
        ufl.as_vector((-(1 / (r * mu)) * A.dx(1), (1 / (r * mu)) * A.dx(0))),
        W.element.interpolation_points()))
    B.name, H.name = "B", "H"

    return MagnetoSolution(mesh, cell_tags, facet_tags, current, r, dx,
                           A, mu, J, B, H)


# ----------------------------------------------------------------------------
# Grandezas globais
# ----------------------------------------------------------------------------
def magnetic_energy(sol):
    """Energia magnética W = pi * ∫ (r/mu) |B|^2 dr dz   (J)."""
    integrando = sol.r / sol.mu * ufl.dot(sol.B, sol.B) * sol.dx
    return np.pi * fem.assemble_scalar(fem.form(integrando))


def inductance(sol):
    """Indutância L = 2 W / I^2   (H)."""
    return 2 * magnetic_energy(sol) / sol.current ** 2


def ohmic_power(sol):
    """Potência dissipada na bobina P = 2 pi ∫ r J^2 / sigma dr dz   (W)."""
    integrando = sol.r * sol.J ** 2 / cfg.fio.sigma * sol.dx(cfg.TAG_COIL)
    return 2 * np.pi * fem.assemble_scalar(fem.form(integrando))


def b_magnitude(sol):
    """|B| por elemento (T)."""
    return np.linalg.norm(sol.B.x.array.reshape(-1, 2), axis=1)


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Teste rápido do solver")
    ap.add_argument("--x-gap-mm", type=float, default=cfg.X_GAP * 1e3)
    ap.add_argument("--current", type=float, default=1.0)
    args = ap.parse_args()
    x_gap = args.x_gap_mm * 1e-3

    mesh, cell_tags, facet_tags = geometry.make_mesh(x_gap)
    sol = solve_magnetostatic(mesh, cell_tags, facet_tags, current=args.current)

    Bmag = b_magnitude(sol)
    B_core = Bmag[cell_tags.find(cfg.TAG_CORE)]
    B_sph = Bmag[cell_tags.find(cfg.TAG_SPH)]
    L = inductance(sol)
    P = ohmic_power(sol)

    print(f"x = {x_gap * 1e3:.3f} mm | I = {args.current} A | "
          f"{mesh.topology.index_map(2).size_local} células")
    print(f"  max |A'|         = {np.abs(sol.A.x.array).max():.3e} Wb")
    print(f"  |B| max no núcleo = {B_core.max():.4f} T")
    print(f"  |B| max na esfera = {B_sph.max():.4f} T")
    print(f"  L                 = {L:.4e} H")
    print(f"  R = P / I^2       = {P / args.current ** 2:.4f} ohm   "
          f"(P = {P:.3f} W)")
    print()
    print("Conferência (valores aproximados, só para checar a ordem de grandeza):")
    print("  I = 1 A, x = 5 mm: |B| max no núcleo ~ 0.08 T, L ~ 2.5e-3 H, R ~ 0.25 ohm")
