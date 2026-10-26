"""
geometry.py
Geometria e malha do levitador axissimétrico usando a API Python do Gmsh.

Plano (r, z): r é o eixo x e z é o eixo y. O eixo de simetria é r = 0.
Só o lado direito (r >= 0) é desenhado.

Regiões:
    1 núcleo, 2 bobina, 3 esfera (semicírculo), 4 ar
Curva marcada:
    10 arco da superfície externa da esfera (usado depois para a força)

Uso rápido (testa só o Gmsh, não precisa do FEniCS):
    python geometry.py                  # checa áreas e gera a malha
    python geometry.py --x-gap-mm 0.5   # outro valor de x
    python geometry.py --gui            # abre a interface gráfica do Gmsh
    python geometry.py --save malha.msh # salva a malha
    python geometry.py --dolfinx        # também converte para o FEniCS
"""
import argparse

import numpy as np
import gmsh

import config as cfg

GDIM = 2


# ----------------------------------------------------------------------------
# 1) Desenho (CAD)
# ----------------------------------------------------------------------------
def _desenhar(x_gap):
    """Cria as 4 formas, funde (fragment) e cria os physical groups."""
    occ = gmsh.model.occ

    # Caixa de ar (grande). As outras formas ficam dentro dela.
    ar = occ.addRectangle(0, cfg.Z_AIR_MIN, 0,
                          cfg.R_AIR, cfg.Z_AIR_MAX - cfg.Z_AIR_MIN)

    # Núcleo: base em z = 0, topo em z = Z_CORE
    nucleo = occ.addRectangle(0, 0, 0, cfg.R_CORE, cfg.Z_CORE)

    # Bobina: ao lado do núcleo (r de R_CORE até R_COIL)
    bobina = occ.addRectangle(cfg.R_CORE, cfg.COIL_Z0, 0,
                              cfg.COIL_THICKNESS, cfg.Z_COIL)

    # Esfera: semicírculo no plano (r,z). O topo fica em z = -x_gap.
    zc = -x_gap - cfg.R_SPH
    p_top = occ.addPoint(0, zc + cfg.R_SPH, 0)
    p_cen = occ.addPoint(0, zc, 0)
    p_bot = occ.addPoint(0, zc - cfg.R_SPH, 0)
    p_dir = occ.addPoint(cfg.R_SPH, zc, 0)
    arco_sup = occ.addCircleArc(p_top, p_cen, p_dir)   # topo -> equador
    arco_inf = occ.addCircleArc(p_dir, p_cen, p_bot)   # equador -> base
    eixo = occ.addLine(p_bot, p_top)                   # fecha no eixo r = 0
    loop = occ.addCurveLoop([arco_sup, arco_inf, eixo])
    esfera = occ.addPlaneSurface([loop])

    # fragment: corta o ar onde estão as outras formas e faz as bordas em
    # comum virarem a MESMA linha (malha conforme entre materiais).
    # out_map[i] diz em quais pedaços cada forma de entrada foi dividida,
    # na ordem em que foram passadas: [ar, nucleo, bobina, esfera].
    ferramentas = [(GDIM, nucleo), (GDIM, bobina), (GDIM, esfera)]
    _, out_map = occ.fragment([(GDIM, ar)], ferramentas)
    occ.synchronize()

    def tags(m):
        return [t for d, t in m if d == GDIM]

    t_nucleo, t_bobina, t_esfera = tags(out_map[1]), tags(out_map[2]), tags(out_map[3])
    usados = set(t_nucleo + t_bobina + t_esfera)
    # o ar é o que sobrou da caixa de ar depois de tirar as outras regiões
    t_ar = [t for t in tags(out_map[0]) if t not in usados]

    # Curvas do arco da esfera: fronteira da esfera que não está no eixo
    fronteira = gmsh.model.getBoundary([(GDIM, t) for t in t_esfera],
                                       oriented=False, combined=False)
    t_arcos = [t for d, t in fronteira
               if occ.getCenterOfMass(d, t)[0] > 1e-9]

    # Physical groups: são os números que viajam junto com a malha pro FEniCS
    pg = gmsh.model.addPhysicalGroup
    pg(GDIM, t_nucleo, tag=cfg.TAG_CORE, name="core")
    pg(GDIM, t_bobina, tag=cfg.TAG_COIL, name="coil")
    pg(GDIM, t_esfera, tag=cfg.TAG_SPH, name="sphere")
    pg(GDIM, t_ar, tag=cfg.TAG_AIR, name="air")
    pg(1, t_arcos, tag=cfg.TAG_SPH_SURF, name="sphere_surface")

    return dict(core=t_nucleo, coil=t_bobina, sphere=t_esfera,
                air=t_ar, sphere_arcs=t_arcos)


# ----------------------------------------------------------------------------
# 2) Tamanho dos elementos (refino no entreferro)
# ----------------------------------------------------------------------------
def _tamanho_malha(x_gap):
    """Elementos pequenos no entreferro, médios no dispositivo, grandes longe."""
    f = gmsh.model.mesh.field

    # tamanho no entreferro: GAP_ELEMS elementos ao longo de x_gap
    h_gap = max(cfg.H_MIN, min(cfg.H_DEV, x_gap / cfg.GAP_ELEMS))

    z_base_esfera = -x_gap - 2 * cfg.R_SPH

    # Caixa 1: região do dispositivo inteiro
    dev = f.add("Box")
    f.setNumber(dev, "VIn", cfg.H_DEV)
    f.setNumber(dev, "VOut", cfg.H_FAR)
    f.setNumber(dev, "XMin", 0.0)
    f.setNumber(dev, "XMax", 1.5 * cfg.R_COIL)
    f.setNumber(dev, "YMin", z_base_esfera - 0.5 * cfg.R_SPH)
    f.setNumber(dev, "YMax", cfg.Z_CORE + 0.5 * cfg.R_CORE)
    f.setNumber(dev, "ZMin", -1.0)
    f.setNumber(dev, "ZMax", 1.0)
    f.setNumber(dev, "Thickness", 2.0 * cfg.R_COIL)

    # Caixa 2: entreferro e arredores (mais fino)
    gap = f.add("Box")
    f.setNumber(gap, "VIn", h_gap)
    f.setNumber(gap, "VOut", cfg.H_FAR)
    f.setNumber(gap, "XMin", 0.0)
    f.setNumber(gap, "XMax", 1.3 * max(cfg.R_CORE, cfg.R_SPH))
    f.setNumber(gap, "YMin", -x_gap - 0.5 * cfg.R_SPH)
    f.setNumber(gap, "YMax", 0.5 * cfg.R_CORE)
    f.setNumber(gap, "ZMin", -1.0)
    f.setNumber(gap, "ZMax", 1.0)
    f.setNumber(gap, "Thickness", 4.0 * h_gap + cfg.R_CORE)

    # vale o menor dos dois campos
    menor = f.add("Min")
    f.setNumbers(menor, "FieldsList", [dev, gap])
    f.setAsBackgroundMesh(menor)

    # só o campo acima manda no tamanho (desliga os outros critérios do Gmsh)
    opt = gmsh.option.setNumber
    opt("Mesh.MeshSizeFromPoints", 0)
    opt("Mesh.MeshSizeFromCurvature", 0)
    opt("Mesh.MeshSizeExtendFromBoundary", 0)
    opt("Mesh.MeshSizeMin", cfg.H_MIN)
    opt("Mesh.MeshSizeMax", cfg.H_FAR)
    opt("Mesh.Algorithm", 6)   # Frontal-Delaunay


# ----------------------------------------------------------------------------
# 3) Montagem completa (assume gmsh.initialize() já chamado)
# ----------------------------------------------------------------------------
def make_gmsh_model(x_gap):
    gmsh.model.add("levitador")
    info = _desenhar(x_gap)
    _tamanho_malha(x_gap)
    gmsh.model.mesh.generate(GDIM)
    return info


def check_geometry(info):
    """Confere áreas e comprimento do arco contra os valores esperados."""
    occ = gmsh.model.occ
    area = lambda ts: sum(occ.getMass(GDIM, t) for t in ts)
    comp = lambda ts: sum(occ.getMass(1, t) for t in ts)

    esperado = {
        "núcleo": (area(info["core"]), cfg.R_CORE * cfg.Z_CORE),
        "bobina": (area(info["coil"]), cfg.COIL_THICKNESS * cfg.Z_COIL),
        "esfera": (area(info["sphere"]), np.pi * cfg.R_SPH ** 2 / 2),
        "arco da esfera": (comp(info["sphere_arcs"]), np.pi * cfg.R_SPH),
    }
    caixa = cfg.R_AIR * (cfg.Z_AIR_MAX - cfg.Z_AIR_MIN)
    outras = sum(v[1] for k, v in esperado.items() if k in ("núcleo", "bobina", "esfera"))
    esperado["ar"] = (area(info["air"]), caixa - outras)

    ok = True
    print(f"{'região':<16}{'gmsh':>14}{'esperado':>14}{'erro rel.':>12}")
    for nome, (obtido, ref) in esperado.items():
        erro = abs(obtido - ref) / ref
        ok &= erro < 1e-6
        print(f"{nome:<16}{obtido:>14.6e}{ref:>14.6e}{erro:>12.1e}")
    return ok


def contar_elementos():
    """Número de elementos 2D por physical group."""
    for _, tag in gmsh.model.getPhysicalGroups(GDIM):
        nome = gmsh.model.getPhysicalName(GDIM, tag)
        n = 0
        for ent in gmsh.model.getEntitiesForPhysicalGroup(GDIM, tag):
            _, _, nodes = gmsh.model.mesh.getElements(GDIM, ent)
            n += sum(len(x) // 3 for x in nodes)   # triângulos: 3 nós cada
        print(f"  {nome:<8} (tag {tag}): {n} triângulos")


# ----------------------------------------------------------------------------
# 4) Interface com o FEniCS
# ----------------------------------------------------------------------------
def make_mesh(x_gap, verbose=False):
    """
    Gera a malha para o valor de x_gap e devolve (mesh, cell_tags, facet_tags)
    no formato do dolfinx. Código pensado para rodar com 1 processo MPI.
    """
    from mpi4py import MPI
    from dolfinx.io import gmshio   # dolfinx 0.9

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
    try:
        make_gmsh_model(x_gap)
        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, MPI.COMM_WORLD, 0, gdim=GDIM)
    finally:
        gmsh.finalize()
    return mesh, cell_tags, facet_tags


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Testa a geometria/malha do TC5")
    ap.add_argument("--x-gap-mm", type=float, default=cfg.X_GAP * 1e3)
    ap.add_argument("--gui", action="store_true", help="abre a janela do Gmsh")
    ap.add_argument("--save", metavar="ARQ.msh", help="salva a malha")
    ap.add_argument("--dolfinx", action="store_true", help="converte pro FEniCS")
    args = ap.parse_args()
    x_gap = args.x_gap_mm * 1e-3

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    info = make_gmsh_model(x_gap)

    print(f"x_gap = {x_gap * 1e3:.3f} mm")
    print("tags após o fragment:", info)
    ok = check_geometry(info)
    print("GEOMETRIA OK" if ok else "ATENÇÃO: alguma área/comprimento não bate")
    print("Malha:")
    contar_elementos()

    if args.save:
        gmsh.write(args.save)
        print("malha salva em", args.save)
    if args.gui:
        gmsh.fltk.run()
    gmsh.finalize()

    if args.dolfinx:
        mesh, cell_tags, facet_tags = make_mesh(x_gap)
        print("dolfinx: células por tag:",
              dict(zip(*np.unique(cell_tags.values, return_counts=True))))
        print("dolfinx: facets do arco da esfera (tag 10):",
              len(facet_tags.find(cfg.TAG_SPH_SURF)))