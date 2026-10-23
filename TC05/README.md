## Organização proposta dos arquivos:

config.py        parâmetros, constantes, classe AWG
geometry.py      gmsh: geometria, refino no gap, physical groups
solver.py        forma fraca, J, μ, solve, B e H
postprocess.py   L (energia), Fz manual, Fz com ds/FacetNormal, B max/méd
parte1_base.py
parte2_parametros.py
parte3_naolinear.py

## Etapas propostas:

Ambiente e smoke test.
Reproduzir o cilindro do prof em .py e conferir contra o valor do FEMM (F ≈ 2,03e-3 N). Isso valida o pipeline antes de mexer na geometria.
Trocar cilindro por esfera, com checks: comprimento do arco = π·r_sph e estudo de convergência de malha.
Parte 1: I = 1 A, campos, Hr e Hz, Fz manual vs ds2, I_eq, fio AWG 18 e potência.
Parte 2: loop dos pontos (mais x_eq), tabelas em CSV, curvas, k.
Parte 3: B×H, pontos de operação, μr (e Picard se der tempo).
Relatório em PDF respondendo aos comentários do prof, mais o .zip.

## Ambiente com conda

Docker não é necessário (fenics), é possível utilizar o conda-forge:

conda create -n fenics-tc5 -c conda-forge fenics-dolfinx=0.9.0 mpich pyvista python-gmsh matplotlib scipy pandas
conda activate fenics-tc5
