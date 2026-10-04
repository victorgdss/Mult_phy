from pathlib import Path
import warnings
import gmsh
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpi4py import MPI
from petsc4py import PETSc
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from dolfinx import fem, geometry
from dolfinx.fem import functionspace, dirichletbc, locate_dofs_topological
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmsh as gmshio
from dolfinx.plot import vtk_mesh
from ufl import TrialFunction, TestFunction, SpatialCoordinate, grad, dot, dx, inner, sqrt
try:
    import pyvista
except ImportError:
    pyvista = None
 
 
# =============================================================================
# 1. PARAMETROS DO ENUNCIADO
# =============================================================================
 
# Geometria - Tabela 1 do TC04
R = 0.12e-3 / 2.0          # raio do tubo capilar [m]
d = 6.0e-3                 # distancia entre a ponta do capilar e o eletrodo [m]
Ra = 2.0e-3                # raio interno do eletrodo [m]
Rb = 9.0e-3                # raio externo do eletrodo [m]
Rc = R                     # raio de curvatura do propelente [m]
 
# Formamida - TC04
gamma = 0.05               # tensao superficial [N/m]
epsilon_r = 84.0           # (referencia; NAO usado no problema, ver nota na secao 4)
epsilon_0 = 8.8541878128e-12  # permissividade do vacuo [F/m]
 
# Goticula - Tabela 2 do TC04
Rg = 0.002e-3              # raio da goticula [m]
rho_formamida = 1.13e3     # densidade [kg/m^3]
 
# Tensao de referencia. Como Laplace e linear, resolve-se uma vez para 1 V.
Vc_ref = 1.0
 
# Parametros auxiliares de geometria/malha (apenas fecham o dominio numerico).
# NOTAS DE MODELAGEM (declarar no relatorio):
#  - z = 0 e fronteira de Neumann (plano de simetria/espelho): a agulha e,
#    efetivamente, mais longa que os 3 mm modelados.
#  - A espessura do eletrodo (2R) e uma escolha numerica, nao do enunciado.
z_nozzle = d / 2.0
z_tip = z_nozzle + R
z_electrode = z_tip + d    # distancia ponta-eletrodo exatamente d
t_electrode = 2.0 * R
z_max = z_electrode + d / 2.0
rho_max = Rb
 
lc_far = 0.25e-3
lc_tip = R / 32.0          # malha principal (P2 converge bem aqui)
 
# Estudo de convergencia de malha (divisores de R para lc_tip)
RODAR_CONVERGENCIA = True
DIVISORES_CONV = (8, 16, 32)       # R/64 falha na localizacao do ponto da ponta no DOLFINx
 
# Tags fisicas do Gmsh
AIR = 1
EMISSOR_BC = 11
ELETRODO_BC = 12
AXIS_BC = 13
OUTER_BC = 14
 
OUTDIR = Path("resultados_tc04")
 
 
# =============================================================================
# 2. FUNCOES AUXILIARES
# =============================================================================
 
def criar_linear_problem(a, L, bcs, prefixo):
    opcoes = {"ksp_type": "preonly", "pc_type": "lu"}
    try:
        return LinearProblem(a, L, bcs=bcs, petsc_options_prefix=prefixo, petsc_options=opcoes)
    except TypeError:
        return LinearProblem(a, L, bcs=bcs, petsc_options=opcoes)
 
 
def converter_modelo_gmsh(comm):
    dados = gmshio.model_to_mesh(gmsh.model, comm, 0, gdim=2)
    if hasattr(dados, "mesh"):
        return dados.mesh, dados.cell_tags, dados.facet_tags
    return dados
 
 
def projetar_escalar(expr, V, prefixo):
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u, v) * dx
    L = inner(expr, v) * dx
    problema = criar_linear_problem(a, L, [], prefixo)
    return problema.solve()
 
 
def _localizar_celulas(msh, pontos):
    arvore = geometry.bb_tree(msh, msh.topology.dim)
    candidatos = geometry.compute_collisions_points(arvore, pontos)
    colididas = geometry.compute_colliding_cells(msh, candidatos, pontos)
    celulas = np.full(len(pontos), -1, dtype=np.int32)
    for i in range(len(pontos)):
        links = colididas.links(i)
        if len(links) > 0:
            celulas[i] = links[0]
    return celulas
 
 
def avaliar_funcao_em_pontos(funcao, pontos_2d):
    msh = funcao.function_space.mesh
    base = np.zeros((len(pontos_2d), 3), dtype=np.float64)
    base[:, :2] = np.asarray(pontos_2d, dtype=np.float64)
 
    pontos = base.copy()
    celulas = _localizar_celulas(msh, pontos)
 
    tentativas = [(dr, dz) for dr in (1.0e-9, 1.0e-8, 1.0e-7) for dz in (0.0, 1.0e-9, 1.0e-8)]
    for dr, dz in tentativas:
        falhos = np.where(celulas < 0)[0]
        if len(falhos) == 0:
            break
        p_try = base[falhos].copy()
        p_try[:, 0] += dr
        p_try[:, 1] += dz
        c_try = _localizar_celulas(msh, p_try)
        ok = c_try >= 0
        pontos[falhos[ok]] = p_try[ok]
        celulas[falhos[ok]] = c_try[ok]
 
    valido = celulas >= 0
    valores = np.full(len(base), np.nan, dtype=float)
    if np.any(valido):
        v = funcao.eval(pontos[valido], celulas[valido])
        valores[valido] = np.asarray(v).reshape(-1).real
    return valores
 
 
def salvar_figura_linha(x, y, xlabel, ylabel, titulo, nome):
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, linewidth=1.8)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(titulo)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / nome, dpi=220)
    plt.close()
 
 
def configurar_pyvista():
    if pyvista is None:
        return False
    pyvista.OFF_SCREEN = True
    try:
        pyvista.start_xvfb(wait=0.1)
    except Exception:
        pass
    return True
 
 
# =============================================================================
# 3. GEOMETRIA E MALHA AXISSIMETRICA
# =============================================================================
 
def construir_malha(comm, lc_tip_loc=None):
    if lc_tip_loc is None:
        lc_tip_loc = lc_tip
    if gmsh.isInitialized():
        gmsh.finalize()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("TC04_eletropropulsor")
 
    # Emissor: parede vertical + ponta circular (calota hemisferica) de raio R.
    p1 = gmsh.model.geo.addPoint(R, 0.0, 0.0, lc_tip_loc)
    p2 = gmsh.model.geo.addPoint(R, z_nozzle, 0.0, lc_tip_loc)
    p3 = gmsh.model.geo.addPoint(0.0, z_tip, 0.0, lc_tip_loc)
    pc = gmsh.model.geo.addPoint(0.0, z_nozzle, 0.0, lc_tip_loc)
 
    # Fronteira externa.
    p4 = gmsh.model.geo.addPoint(rho_max, 0.0, 0.0, lc_far)
    p5 = gmsh.model.geo.addPoint(rho_max, z_max, 0.0, lc_far)
    p6 = gmsh.model.geo.addPoint(0.0, z_max, 0.0, lc_far)
 
    # Eletrodo anular.
    p7 = gmsh.model.geo.addPoint(Ra, z_electrode, 0.0, lc_far)
    p8 = gmsh.model.geo.addPoint(rho_max, z_electrode, 0.0, lc_far)
    p9 = gmsh.model.geo.addPoint(rho_max, z_electrode + t_electrode, 0.0, lc_far)
    p10 = gmsh.model.geo.addPoint(Ra, z_electrode + t_electrode, 0.0, lc_far)
 
    l_em_wall = gmsh.model.geo.addLine(p1, p2)
    l_em_arc = gmsh.model.geo.addCircleArc(p2, pc, p3)
 
    l_el_bottom = gmsh.model.geo.addLine(p7, p8)
    l_el_top = gmsh.model.geo.addLine(p9, p10)
    l_el_left = gmsh.model.geo.addLine(p10, p7)
 
    l_air_bottom = gmsh.model.geo.addLine(p1, p4)
    l_air_right_low = gmsh.model.geo.addLine(p4, p8)
    l_air_right_high = gmsh.model.geo.addLine(p9, p5)
    l_air_top = gmsh.model.geo.addLine(p5, p6)
    l_air_axis = gmsh.model.geo.addLine(p6, p3)
 
    air_loop = gmsh.model.geo.addCurveLoop([
        l_air_bottom, l_air_right_low,
        -l_el_bottom, -l_el_left, -l_el_top,
        l_air_right_high, l_air_top, l_air_axis,
        -l_em_arc, -l_em_wall,
    ])
    air_surface = gmsh.model.geo.addPlaneSurface([air_loop])
    gmsh.model.geo.synchronize()
 
    gmsh.model.addPhysicalGroup(2, [air_surface], AIR)
    gmsh.model.setPhysicalName(2, AIR, "dominio")
    gmsh.model.addPhysicalGroup(1, [l_em_wall, l_em_arc], EMISSOR_BC)
    gmsh.model.setPhysicalName(1, EMISSOR_BC, "emissor")
    gmsh.model.addPhysicalGroup(1, [l_el_bottom, l_el_left, l_el_top], ELETRODO_BC)
    gmsh.model.setPhysicalName(1, ELETRODO_BC, "eletrodo")
    gmsh.model.addPhysicalGroup(1, [l_air_axis], AXIS_BC)
    gmsh.model.setPhysicalName(1, AXIS_BC, "eixo")
    gmsh.model.addPhysicalGroup(1, [l_air_bottom, l_air_right_low, l_air_right_high, l_air_top], OUTER_BC)
    gmsh.model.setPhysicalName(1, OUTER_BC, "fronteira_externa")
 
    # Refinamento na ponta do emissor e na borda interna do eletrodo.
    campo_dist_em = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(campo_dist_em, "CurvesList", [l_em_arc])
    campo_thr_em = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(campo_thr_em, "InField", campo_dist_em)
    gmsh.model.mesh.field.setNumber(campo_thr_em, "SizeMin", lc_tip_loc)
    gmsh.model.mesh.field.setNumber(campo_thr_em, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(campo_thr_em, "DistMin", 2.0 * R)
    gmsh.model.mesh.field.setNumber(campo_thr_em, "DistMax", 1.0e-3)
 
    campo_dist_el = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(campo_dist_el, "CurvesList", [l_el_left])
    campo_thr_el = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(campo_thr_el, "InField", campo_dist_el)
    gmsh.model.mesh.field.setNumber(campo_thr_el, "SizeMin", 0.05e-3)
    gmsh.model.mesh.field.setNumber(campo_thr_el, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(campo_thr_el, "DistMin", 0.2e-3)
    gmsh.model.mesh.field.setNumber(campo_thr_el, "DistMax", 1.0e-3)
 
    campo_min = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(campo_min, "FieldsList", [campo_thr_em, campo_thr_el])
    gmsh.model.mesh.field.setAsBackgroundMesh(campo_min)
 
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_far)
    gmsh.model.mesh.generate(2)
 
    try:
        msh, cell_tags, facet_tags = converter_modelo_gmsh(comm)
    finally:
        gmsh.finalize()
 
    return msh, cell_tags, facet_tags
 
 
# =============================================================================
# 4. PROBLEMA ELETROSTATICO FEM
# =============================================================================
 
def resolver_eletrostatica(msh, facet_tags, prefixo="tc04"):
    Vh = functionspace(msh, ("Lagrange", 2))
    fdim = msh.topology.dim - 1
 
    facetas_emissor = facet_tags.find(EMISSOR_BC)
    facetas_eletrodo = facet_tags.find(ELETRODO_BC)
    dofs_emissor = locate_dofs_topological(Vh, fdim, facetas_emissor)
    dofs_eletrodo = locate_dofs_topological(Vh, fdim, facetas_eletrodo)
 
    bc_emissor = dirichletbc(PETSc.ScalarType(+Vc_ref / 2.0), dofs_emissor, Vh)
    bc_eletrodo = dirichletbc(PETSc.ScalarType(-Vc_ref / 2.0), dofs_eletrodo, Vh)
 
    phi = TrialFunction(Vh)
    w = TestFunction(Vh)
 
    x = SpatialCoordinate(msh)
    r = x[0]
    epsilon = fem.Constant(msh, PETSc.ScalarType(1.0))
 
    a = epsilon * r * dot(grad(w), grad(phi)) * dx
    L = fem.Constant(msh, PETSc.ScalarType(0.0)) * w * dx
 
    problema = criar_linear_problem(a, L, [bc_emissor, bc_eletrodo], f"{prefixo}_phi_")
    phi_h = problema.solve()
    phi_h.name = "phi"
 
    # Com phi em P2, grad(phi) e P1 por celula: a projecao em DG1 e local e
    # exata (sem suavizacao entre elementos, que subestimava E na ponta).
    Vdg1 = functionspace(msh, ("DG", 1))
    Ez_h = projetar_escalar(-grad(phi_h)[1], Vdg1, f"{prefixo}_ez_")
    Er_h = projetar_escalar(-grad(phi_h)[0], Vdg1, f"{prefixo}_er_")
    Ez_h.name = "Ez"
    Er_h.name = "Er"
 
    # Modulo de E para visualizacao (DG0).
    Vdg0 = functionspace(msh, ("DG", 0))
    Emag_h = projetar_escalar(sqrt(inner(grad(phi_h), grad(phi_h))), Vdg0, f"{prefixo}_emag_")
    Emag_h.name = "E_mag"
 
    return phi_h, Ez_h, Er_h, Emag_h
 
 
# =============================================================================
# 5. POS-PROCESSAMENTO NO EIXO E TENSAO CRITICA
# =============================================================================
 
def perfil_no_eixo(Ez_h, n=700):
    # Avalia SOMENTE no eixo r = 0. Amostragem logaritmica: o campo varia na
    # escala de R perto da ponta.
    ds = max(1.0e-9, R * 1.0e-4)
    s = np.concatenate(([ds], np.geomspace(1.0e-6, d - ds, n - 1)))
    z = z_tip + s
    pontos = np.column_stack((np.zeros_like(z), z))
    Ez = avaliar_funcao_em_pontos(Ez_h, pontos)
 
    if not np.isfinite(Ez[0]):
        raise RuntimeError("Nao foi possivel avaliar o campo na ponta do emissor (primeiro ponto do eixo).")
    if np.any(~np.isfinite(Ez)):
        bons = np.isfinite(Ez)
        warnings.warn(f"{int((~bons).sum())} ponto(s) do eixo nao localizado(s); interpolados dos vizinhos.")
        Ez = np.interp(s, s[bons], Ez[bons])
 
    # Emissor em maior potencial: Ez deve apontar para +z (eletrodo).
    assert np.median(Ez) > 0.0, "Ez com sinal inesperado: revisar condicoes de contorno."
 
    return s, Ez, ds
 
 
def calcular_tensoes_criticas(E_tip_1V):
    E_crit = np.sqrt(np.pi * gamma / (2.0 * epsilon_0 * R))                       # Eq. (1)
    V_crit_analitico = np.sqrt(gamma * Rc / epsilon_0) * np.log(4.0 * d / Rc)     # Eq. (2)
    V_crit_sim = E_crit / E_tip_1V                                                # linearidade
    erro = abs(V_crit_sim - V_crit_analitico) / V_crit_analitico * 100.0
    return E_crit, V_crit_sim, V_crit_analitico, erro
 
 
def verificar_campo_na_calota(Ez_h, Er_h, V_crit, E_crit, n=90):
    eps = 1.0e-8
    th = np.linspace(0.0, np.pi / 2.0, n)
    r_p = (R + eps) * np.sin(th)
    z_p = z_nozzle + (R + eps) * np.cos(th)
    pts = np.column_stack((r_p, z_p))
    Ez = avaliar_funcao_em_pontos(Ez_h, pts) * V_crit
    Er = avaliar_funcao_em_pontos(Er_h, pts) * V_crit
    En = Er * np.sin(th) + Ez * np.cos(th)
 
    plt.figure(figsize=(8, 5))
    plt.plot(np.degrees(th), En, linewidth=1.8, label="E normal na calota (V = V_crit simulada)")
    plt.axhline(E_crit, linestyle="--", linewidth=1.2, label="Campo critico (Eq. 1)")
    plt.xlabel("Angulo a partir do eixo (graus)  [0 = apice, 90 = ombro]")
    plt.ylabel("Campo eletrico normal (V/m)")
    plt.title("Verificacao da condicao de campo critico na calota")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTDIR / "08_verificacao_campo_critico.png", dpi=220)
    plt.close()
    # Maximo robusto: filtro de mediana (janela 5) remove picos isolados causados
    # pela descontinuidade de DG1 entre elementos vizinhos da calota.
    pad = np.pad(En, 2, mode="edge")
    En_filtrado = np.median(np.lib.stride_tricks.sliding_window_view(pad, 5), axis=1)
    return float(np.nanmax(En_filtrado)), float(En[0])
 
 
def estudo_convergencia(comm):
    E_crit = np.sqrt(np.pi * gamma / (2.0 * epsilon_0 * R))
    linhas = []
    print("\n--- Convergencia de malha (P2) ---")
    for div in DIVISORES_CONV:
        msh_i, _, ft_i = construir_malha(comm, R / div)
        _, Ez_i, _, _ = resolver_eletrostatica(msh_i, ft_i, prefixo=f"conv{div}")
        try:
            _, Ez_prof, _ = perfil_no_eixo(Ez_i, n=50)
        except RuntimeError as err:
            print(f"lc_tip = R/{div:<3d}  falhou: {err}")
            continue
        Et = abs(Ez_prof[0])
        ncel = msh_i.topology.index_map(msh_i.topology.dim).size_local
        Vc = E_crit / Et
        print(f"lc_tip = R/{div:<3d}  celulas = {ncel:<8d}  E_tip(1V) = {Et:10.1f} V/m   V_crit = {Vc:8.2f} V")
        linhas.append((div, ncel, Et, Vc))
    np.savetxt(OUTDIR / "convergencia_malha.csv", np.array(linhas), delimiter=",",
               header="div_R,celulas,E_tip_1V_Vpm,V_crit_V", comments="")
    return linhas
 
 
# =============================================================================
# 6. MODELO DINAMICO DA GOTICULA
# =============================================================================
 
def resolver_dinamica(s, Ez_1V, V_crit_sim):
    Q = 8.0 * np.pi * np.sqrt(epsilon_0 * gamma) * Rg ** 1.5     # Eq. (8)
    m = (4.0 / 3.0) * np.pi * Rg ** 3 * rho_formamida
 
    Ez_crit = Ez_1V * V_crit_sim
    interp_Ez = PchipInterpolator(s, Ez_crit, extrapolate=True)
 
    # O centro da goticula parte a uma distancia Rg da ponta (ela tem raio Rg).
    s_ini = float(max(Rg, s[0]))
    s_fim = float(s[-1])
 
    def dinamica(t, y):
        pos = float(np.clip(y[0], s_ini, s_fim))
        vel = y[1]
        Fe = Q * float(interp_Ez(pos))
        return [vel, Fe / m]          # Eqs. (9)-(10)
 
    def chegou_ao_eletrodo(t, y):
        return y[0] - s_fim
 
    chegou_ao_eletrodo.terminal = True
    chegou_ao_eletrodo.direction = 1
 
    y0 = [s_ini, 0.0]   # velocidade inicial nula; aceleracao inicial vem da Eq. (10)
 
    sol = solve_ivp(dinamica, (0.0, 1.0e-3), y0, method="RK45", events=chegou_ao_eletrodo,
                    max_step=1.0e-7, rtol=1.0e-8, atol=1.0e-11)
 
    if not sol.success:
        raise RuntimeError(f"Falha na integracao da goticula: {sol.message}")
    if len(sol.t_events[0]) == 0:
        raise RuntimeError("A goticula nao atingiu o eletrodo no intervalo de integracao.")
 
    return Q, m, Ez_crit, interp_Ez, sol, s_ini, s_fim
 
 
# =============================================================================
# 7. FIGURAS
# =============================================================================
 
def plotar_campos(msh, phi_h, Emag_h):
    if not configurar_pyvista():
        warnings.warn("PyVista nao esta instalado; figuras 2D serao ignoradas.")
        return
 
    cam = "viridis"
 
    # Malha completa.
    grid_mesh = pyvista.UnstructuredGrid(*vtk_mesh(msh, msh.topology.dim))
    p = pyvista.Plotter(off_screen=True, window_size=(1200, 800))
    p.add_mesh(grid_mesh, show_edges=True)
    pontos_rotulo = np.array([[0.0, z_tip, 0.0], [Ra, z_electrode, 0.0]])
    p.add_point_labels(pontos_rotulo, ["Ponta do emissor", "Eletrodo"], point_size=5, font_size=15)
    p.view_xy()
    p.screenshot(str(OUTDIR / "01_malha_geometria.png"))
    p.close()
 
    # Zoom da malha na ponta.
    p = pyvista.Plotter(off_screen=True, window_size=(1200, 800))
    p.add_mesh(grid_mesh, show_edges=True)
    p.view_xy()
    p.camera.focal_point = (0.0, z_tip, 0.0)
    p.camera.position = (0.0, z_tip, 1.0)
    p.camera.parallel_projection = True
    p.camera.parallel_scale = 4.0 * R
    p.screenshot(str(OUTDIR / "01b_malha_zoom_ponta.png"))
    p.close()
 
    # Potencial: interpola para P1 (mapeamento direto malha-pontos).
    msh_p1 = fem.Function(functionspace(msh, ("Lagrange", 1)))
    msh_p1.interpolate(phi_h)
    top, tipos, geo = vtk_mesh(msh_p1.function_space)
    grid_phi = pyvista.UnstructuredGrid(top, tipos, geo)
    grid_phi.point_data["phi"] = np.asarray(msh_p1.x.array).real
    p = pyvista.Plotter(off_screen=True, window_size=(1200, 800))
    p.add_mesh(grid_phi, scalars="phi", cmap=cam, show_edges=False,
               scalar_bar_args={"title": "phi [V]", "vertical": True})
    p.add_mesh(grid_phi.contour(isosurfaces=15, scalars="phi"), color="white", line_width=1)
    p.view_xy()
    p.screenshot(str(OUTDIR / "02_potencial_1V.png"))
    p.close()
 
    # Modulo do campo (escala logaritmica).
    top, tipos, geo = vtk_mesh(msh, msh.topology.dim)
    grid_E = pyvista.UnstructuredGrid(top, tipos, geo)
    ncel = msh.topology.index_map(msh.topology.dim).size_local
    grid_E.cell_data["E_mag"] = np.maximum(np.asarray(Emag_h.x.array).real[:ncel], 1.0e-3)
    p = pyvista.Plotter(off_screen=True, window_size=(1200, 800))
    p.add_mesh(grid_E, scalars="E_mag", cmap=cam, log_scale=True, show_edges=False,
               scalar_bar_args={"title": "|E| [V/m] (1 V)", "vertical": True})
    p.view_xy()
    p.screenshot(str(OUTDIR / "03_modulo_campo_1V.png"))
    p.close()
 
 
def plotar_resultados_1d(s, Ez_1V, Q, interp_Ez, sol, s_ini, s_fim):
    s_mm = s * 1.0e3
 
    salvar_figura_linha(s_mm, Ez_1V, "Distancia a partir da ponta do emissor (mm)",
                        "Ez para Vc = 1 V (V/m)", "Campo eletrico axial no eixo de simetria",
                        "04_Ez_eixo_1V.png")
 
    s_dense = np.geomspace(s_ini, s_fim, 1500)
    Fe_dense = Q * interp_Ez(s_dense)
    salvar_figura_linha(s_dense * 1.0e3, Fe_dense, "Distancia a partir da ponta do emissor (mm)",
                        "Forca eletrica Fe (N)", "Forca eletrica sobre a goticula",
                        "05_forca_posicao.png")
 
    salvar_figura_linha(sol.y[0] * 1.0e3, sol.y[1], "Distancia a partir da ponta do emissor (mm)",
                        "Velocidade da goticula (m/s)", "Velocidade da goticula em funcao da posicao",
                        "06_velocidade_posicao.png")
 
    salvar_figura_linha(sol.t * 1.0e6, sol.y[0] * 1.0e3, "Tempo (us)",
                        "Distancia a partir da ponta do emissor (mm)",
                        "Posicao da goticula em funcao do tempo", "07_posicao_tempo.png")
 
 
# =============================================================================
# 8. PROGRAMA PRINCIPAL
# =============================================================================
 
def main():
    comm = MPI.COMM_WORLD
    if comm.size != 1:
        raise RuntimeError("Este script foi preparado para execucao serial: use 'python3 TC04.py'.")
 
    OUTDIR.mkdir(parents=True, exist_ok=True)
 
    print("=" * 72)
    print("TC04 - ELETROPROPULSOR AXISSIMETRICO")
    print("=" * 72)
 
    if RODAR_CONVERGENCIA:
        estudo_convergencia(comm)
 
    msh, cell_tags, facet_tags = construir_malha(comm)
    print(f"\nMalha principal (lc_tip = R/{R / lc_tip:.0f}): "
          f"{msh.topology.index_map(msh.topology.dim).size_local} celulas")
 
    phi_h, Ez_h, Er_h, Emag_h = resolver_eletrostatica(msh, facet_tags)
    print(f"Potencial minimo: {np.min(np.asarray(phi_h.x.array).real):.6g} V")
    print(f"Potencial maximo: {np.max(np.asarray(phi_h.x.array).real):.6g} V")
 
    s, Ez_1V, ds = perfil_no_eixo(Ez_h)
    E_tip_1V = abs(Ez_1V[0])
    E_crit, V_crit_sim, V_crit_analitico, erro = calcular_tensoes_criticas(E_tip_1V)
 
    print("\n--- Campo critico e tensoes ---")
    print(f"Primeiro ponto do eixo: {ds * 1e6:.6f} um apos a ponta")
    print(f"Ez na ponta para Vc = 1 V: {E_tip_1V:.6e} V/m")
    print(f"Campo critico (Eq. 1): {E_crit:.6e} V/m")
    print(f"Tensao critica simulada: {V_crit_sim:.6f} V")
    print(f"Tensao critica analitica (Eq. 2): {V_crit_analitico:.6f} V")
    print(f"Erro relativo: {erro:.4f} %")
 
    En_max, En_apice = verificar_campo_na_calota(Ez_h, Er_h, V_crit_sim, E_crit)
    print(f"Verificacao na calota (V = V_crit, max filtrado): max(En)/E_crit = {En_max / E_crit:.4f}; "
          f"En(apice)/E_crit = {En_apice / E_crit:.4f}")
 
    Q, m, Ez_crit, interp_Ez, sol, s_ini, s_fim = resolver_dinamica(s, Ez_1V, V_crit_sim)
 
    v_saida = sol.y[1, -1]
    v_energia = np.sqrt(2.0 * Q / m * float(interp_Ez.integrate(s_ini, s_fim)))
 
    print("\n--- Goticula ---")
    print(f"Carga Q (Eq. 8): {Q:.6e} C")
    print(f"Massa: {m:.6e} kg")
    print(f"Posicao inicial do centro: {s_ini * 1e6:.3f} um")
    print(f"Tempo de travessia: {sol.t[-1]:.6e} s")
    print(f"Velocidade de saida (RK45): {v_saida:.6f} m/s")
    print(f"Velocidade de saida (balanco de energia, conferencia): {v_energia:.6f} m/s")
 
    Fe_crit = Q * Ez_crit
    dados = np.column_stack((s, Ez_1V, Ez_crit, Fe_crit))
    np.savetxt(OUTDIR / "perfil_eixo.csv", dados, delimiter=",",
               header="s_m,Ez_1V_Vpm,Ez_crit_Vpm,Fe_N", comments="")
 
    with open(OUTDIR / "resumo_resultados.txt", "w", encoding="utf-8") as f:
        f.write("TC04 - resumo dos resultados\n")
        f.write(f"E_crit = {E_crit:.12e} V/m\n")
        f.write(f"E_tip_1V = {E_tip_1V:.12e} V/m\n")
        f.write(f"V_crit_sim = {V_crit_sim:.12e} V\n")
        f.write(f"V_crit_analitico = {V_crit_analitico:.12e} V\n")
        f.write(f"erro_relativo = {erro:.8f} %\n")
        f.write(f"max_En_calota/E_crit = {En_max / E_crit:.8f}\n")
        f.write(f"Q = {Q:.12e} C\n")
        f.write(f"m = {m:.12e} kg\n")
        f.write(f"s_inicial = {s_ini:.12e} m\n")
        f.write(f"tempo_travessia = {sol.t[-1]:.12e} s\n")
        f.write(f"velocidade_saida = {v_saida:.12e} m/s\n")
        f.write(f"velocidade_saida_energia = {v_energia:.12e} m/s\n")
 
    plotar_campos(msh, phi_h, Emag_h)
    plotar_resultados_1d(s, Ez_1V, Q, interp_Ez, sol, s_ini, s_fim)
 
    print(f"\nArquivos salvos em: {OUTDIR.resolve()}")
 
 
if __name__ == "__main__":
    main()