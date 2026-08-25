import os
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
Bl0 = 4.95      # fator de força (N/A) em repouso, usado no modelo linear
L  = 266e-6     # indutância da bobina (H)
R  = 3.3        # resistência da bobina (ohm)

# =========================================================
# 2) LER O ÁUDIO GRAVADO E TRANSFORMAR EM Vin(t)
# =========================================================

pasta_script = os.path.dirname(os.path.abspath(__file__))
pasta_saida = os.path.join(pasta_script, "resultados")
os.makedirs(pasta_saida, exist_ok=True)

caminho_audio = os.path.join(pasta_script, "audio.wav")

fs, audio = wavfile.read(caminho_audio)   # fs = taxa de amostragem (Hz)

# Se o áudio for estéreo (2 canais), pega só um canal
if audio.ndim > 1:
    audio = audio[:, 0]

# Converte para float e normaliza entre -1 e 1
audio = audio.astype(np.float64)
audio = audio / np.max(np.abs(audio))

# Escala para satisfazer Vin < 2.0V
V_PICO = 1.8
Vin_amostras = V_PICO * audio

# Vetor de tempo correspondente a cada amostra
t_audio = np.arange(len(Vin_amostras)) / fs

# Cria uma FUNÇÃO u(t) que o integrador pode chamar em qualquer instante
# (o integrador não usa passos fixos, então precisamos poder avaliar
#  Vin em tempos "quebrados" entre as amostras -> interpolação)

u_interp = interp1d(t_audio, Vin_amostras, bounds_error=False, fill_value=0.0)

# =========================================================
# 3) MODELO NÃO LINEAR DE Bl(x)  (Atividade 2/3)
# =========================================================

# Excursão máxima do cone no modelo linear (você deve estimar isso
# simulando primeiro o modelo LINEAR e vendo o pico de x(t) obtido).
# Aqui deixo um valor de EXEMPLO - troque pelo valor real que você medir!

x_max_linear = 3e-3   # <-- SUBSTITUA pelo pico real de x(t) da simulação linear


x1 = 0.75 * x_max_linear   # até aqui, Bl constante
x2 = 1.50 * x_max_linear   # a partir daqui, Bl = 0

def Bl_de_x(x):
    """
    Fator de força em função da posição do cone.
    |x| <= x1        -> Bl = Bl0 (constante, igual ao modelo linear)
    x1 <= |x| <= x2  -> decaimento polinomial de ordem 2, de Bl0 até 0
    |x| >= x2        -> Bl = 0
    """
    ax = abs(x)
    if ax <= x1:
        return Bl0
    elif ax >= x2:
        return 0.0
    else:
        # polinômio de ordem 2 que vale Bl0 em x1 e 0 em x2,
        # com derivada nula em x2 (transição suave, "encosta" em zero)
        frac = (ax - x1) / (x2 - x1)   # vai de 0 (em x1) a 1 (em x2)
        return Bl0 * (1 - frac)**2

# =========================================================
# 4) EQUAÇÕES DE ESTADO  z = [i, x, v]
# =========================================================

NAO_LINEAR = False   # <-- False para rodar a Atividade 1 (modelo linear), True caso contrário

def derivadas(t, z):
    i, x, v = z
    u = u_interp(t)  # tensão de entrada no instante t

    Bl = Bl_de_x(x) if NAO_LINEAR else Bl0

    di_dt = -(R/L)*i - (Bl/L)*v + u/L
    dx_dt = v
    dv_dt = (Bl/m)*i - (k/m)*x - (b/m)*v

    return [di_dt, dx_dt, dv_dt]

# =========================================================
# 5) RESOLVER A ODE (condições iniciais nulas)
# =========================================================

z0 = [0.0, 0.0, 0.0]     # i(0)=0, x(0)=0, v(0)=0
t_span = (t_audio[0], t_audio[-1])

# t_eval: pontos onde queremos a solução (mesma taxa do áudio, para facilitar salvar depois)
sol = solve_ivp(derivadas, t_span, z0, t_eval=t_audio, method="RK45", max_step=1/fs)

i_t = sol.y[0]
x_t = sol.y[1]
v_t = sol.y[2]

x_max_linear = np.max(np.abs(x_t))
print(f"Excursão máxima do cone (modelo linear): {x_max_linear:.6e} m")

# Aceleração ẍ(t): derivamos a 3ª equação de estado diretamente
# (mais preciso que derivar v_t numericamente)
Bl_vetor = np.array([Bl_de_x(x) if NAO_LINEAR else Bl0 for x in x_t])
a_t = (Bl_vetor/m)*i_t - (k/m)*x_t - (b/m)*v_t   # isso é ẍ(t)

# =========================================================
# 6) GRÁFICOS NO TEMPO
# =========================================================

fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharex=True)
axs[0].plot(t_audio, Vin_amostras); axs[0].set_ylabel("Vin (V)")
axs[1].plot(t_audio, i_t);          axs[1].set_ylabel("i (A)")
axs[2].plot(t_audio, x_t);          axs[2].set_ylabel("x (m)")
axs[3].plot(t_audio, a_t);          axs[3].set_ylabel("ẍ (m/s²)")
axs[3].set_xlabel("tempo (s)")
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "resposta_temporal.png"))
plt.show()

# =========================================================
# 7) ESPECTROS (FFT) DE Vin E ẍ, PARA COMPARAR
# =========================================================

def espectro(sinal, fs):
    N = len(sinal)
    freqs = np.fft.rfftfreq(N, d=1/fs)
    mag = np.abs(np.fft.rfft(sinal)) / N
    return freqs, mag

f_vin, mag_vin = espectro(Vin_amostras, fs)
f_a, mag_a = espectro(a_t, fs)

fig, axs = plt.subplots(2, 1, figsize=(10, 6))
axs[0].plot(f_vin, mag_vin); axs[0].set_title("Espectro de Vin")
axs[0].set_xlim(0, 5000)
axs[1].plot(f_a, mag_a); axs[1].set_title("Espectro de ẍ")
axs[1].set_xlim(0, 5000)
for ax in axs:
    ax.set_xlabel("Frequência (Hz)")
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "espectros.png"))
plt.show()

# =========================================================
# 8) SALVAR ẍ(t) COMO ARQUIVO DE ÁUDIO
# =========================================================

#  Normalizar de volta para o formato de áudio (inteiro 16 bits)

a_norm = a_t / np.max(np.abs(a_t))          # normaliza entre -1 e 1
a_int16 = (a_norm * 32767).astype(np.int16) # converte para int16

wavfile.write(os.path.join(pasta_saida, "TC02-out1.wav"), fs, a_int16)

print("Simulação concluída. Arquivos gerados: resposta_temporal.png, espectros.png, TC02-out1.wav")
