import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# =========================================================
# 1) PARÂMETROS DO ALTO-FALANTE (dados no slide 14)
# =========================================================
m  = 14.35e-3   # massa do cone+bobina (kg)
b  = 0.786      # amortecimento mecânico (kg/s)
k  = 1852       # rigidez da suspensão (N/m)
Bl0 = 4.95      # fator de força EM REPOUSO (N/A) - usado no modelo linear
L  = 266e-6     # indutância da bobina (H)
R  = 3.3        # resistência da bobina (ohm)

# =========================================================
# 2) LER O ÁUDIO GRAVADO E TRANSFORMAR EM Vin(t)
# =========================================================
fs, audio = wavfile.read("TC02-in.wav")   # fs = taxa de amostragem (Hz), audio = amostras

# Se o áudio for estéreo (2 canais), pega só um canal
if audio.ndim > 1:
    audio = audio[:, 0]

# Converte para float e normaliza entre -1 e 1
audio = audio.astype(np.float64)
audio = audio / np.max(np.abs(audio))

# Escala para satisfazer Vin < 2.0V (aqui uso 1.8V de margem de segurança)
V_PICO = 1.8
Vin_amostras = V_PICO * audio

# Vetor de tempo correspondente a cada amostra
t_audio = np.arange(len(Vin_amostras)) / fs

# Cria uma FUNÇÃO u(t) que o integrador pode chamar em qualquer instante
# (o integrador não usa passos fixos, então precisamos poder avaliar
#  Vin em tempos "quebrados" entre as amostras -> interpolação)
u_interp = interp1d(t_audio, Vin_amostras, bounds_error=False, fill_value=0.0)