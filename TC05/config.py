"""
config.py
Parâmetros do TC5 (levitador magnético axissimétrico).

Tudo em unidades SI. Os outros arquivos fazem `import config as cfg`,
então qualquer mudança de parâmetro é feita só aqui.
"""
import numpy as np

# ---------------------------------------------------------------- constantes
MU0 = 4e-7 * np.pi        # permeabilidade do vácuo (H/m)
G = 9.81                  # gravidade (m/s^2)

# -------------------------------------------------------------------- núcleo
R_CORE = 0.75e-2          # raio do núcleo (m)
Z_CORE = 6.0e-2           # altura do núcleo (m)
MUR_CORE = 1200           # permeabilidade relativa (aproximação linear)

# -------------------------------------------------------------------- esfera
R_SPH = 1.0e-2            # raio da esfera (m)
M_SPH = 22.0e-3           # massa (kg)
MUR_SPH = 1200            # permeabilidade relativa (aproximação linear)
X_GAP = 0.5e-2            # distância núcleo-esfera nominal ("x" do enunciado) (m)
PESO = M_SPH * G          # força peso (N)

# ---------------------------------------------------------------------- fio
class AWG:
    """Fio de cobre padrão AWG (mesma fórmula do notebook da Aula 18)."""
    def __init__(self, tipo):
        self.tipo = tipo
        self.diameter = 0.127e-3 * 92 ** ((36 - tipo) / 39)   # (m)
        self.area = np.pi / 4 * self.diameter ** 2            # (m^2)
        self.sigma = 58.0e6                                   # condutividade do cobre (S/m)


AWG_TIPO = 18
fio = AWG(AWG_TIPO)

# ------------------------------------------------------------------- bobina
N_TURNS = 200             # número de espiras
Z_COIL = 4.0e-2           # altura da bobina (m)

# Espessura radial "equivalente": fator de enchimento = 1, ou seja,
# área da bobina no plano (r,z) = N * área do fio
COIL_THICKNESS = N_TURNS * fio.area / Z_COIL
R_COIL = R_CORE + COIL_THICKNESS

# DECISÃO EM ABERTO: o enunciado não diz em que altura do núcleo a bobina fica.
# Aqui: centralizada. Para encostar a bobina na base do núcleo use COIL_Z0 = 0.
# Se mudar, descrever e justificar no relatório.
COIL_Z0 = (Z_CORE - Z_COIL) / 2   # z da base da bobina (m), com z=0 na base do núcleo

# ------------------------------------------------------------- caixa de ar
# Mesmos valores da Aula 18. Depois vale testar caixas maiores (sensibilidade).
R_AIR = 5 * R_COIL
Z_AIR = 2.5 * Z_CORE
Z_AIR_MIN = -Z_AIR / 2
Z_AIR_MAX = Z_AIR / 2

# ----------------------------------------------- tags (physical groups do Gmsh)
TAG_CORE = 1              # superfície: núcleo
TAG_COIL = 2              # superfície: bobina
TAG_SPH = 3               # superfície: esfera
TAG_AIR = 4               # superfície: ar
TAG_SPH_SURF = 10         # curva: arco da superfície externa da esfera

# -------------------------------------------------------------------- malha
H_DEV = R_CORE / 8        # tamanho dos elementos na região do dispositivo (m)
H_FAR = H_DEV             # tamanho dos elementos longe do dispositivo (m)
GAP_ELEMS = 4             # nº de elementos ao longo da altura do entreferro
H_MIN = 0.05e-3           # tamanho mínimo permitido (m)
