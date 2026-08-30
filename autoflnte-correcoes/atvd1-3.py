import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

m = 14.35e-3
b = 0.786
k = 1852
Bl0 = 4.95
L = 266e-6
R = 3.3

pasta_script = os.path.dirname(os.path.abspath(__file__))
pasta_saida = os.path.join(pasta_script, "resultados")
os.makedirs(pasta_saida, exist_ok=True)

caminho_audio = os.path.join(pasta_script, "TC02-in.wav")
fs, audio = wavfile.read(caminho_audio)

if audio.ndim > 1:
    audio = audio[:, 0]

audio = audio.astype(np.float64)
audio = audio / np.max(np.abs(audio))

V_PICO = 1.8
Vin_amostras = V_PICO * audio

t_audio = np.arange(len(Vin_amostras)) / fs

u_interp = interp1d(t_audio, Vin_amostras, bounds_error=False, fill_value=0.0)

x_max_linear = 4.3514e-4

x1 = 0.75 * x_max_linear
x2 = 1.50 * x_max_linear

def Bl_de_x(x):
    ax = abs(x)
    if ax <= x1:
        return Bl0
    elif ax >= x2:
        return 0.0
    else:
        frac = (ax - x1) / (x2 - x1)
        return Bl0 * (1 - frac)**2

NAO_LINEAR = True
label_modelo = "nao_linear" if NAO_LINEAR else "linear"

def derivadas(t, z):
    i, x, v = z
    u = u_interp(t)

    Bl = Bl_de_x(x) if NAO_LINEAR else Bl0

    di_dt = -(R/L)*i - (Bl/L)*v + u/L
    dx_dt = v
    dv_dt = (Bl/m)*i - (k/m)*x - (b/m)*v

    return [di_dt, dx_dt, dv_dt]

z0 = [0.0, 0.0, 0.0]
t_span = (t_audio[0], t_audio[-1])

sol = solve_ivp(derivadas, t_span, z0, t_eval=t_audio, method="RK45", max_step=1/fs)

i_t = sol.y[0]
x_t = sol.y[1]
v_t = sol.y[2]

x_max_linear = np.max(np.abs(x_t))
print(f"Excursão máxima do cone (modelo {label_modelo}): {x_max_linear:.6e} m")

Bl_vetor = np.array([Bl_de_x(x) if NAO_LINEAR else Bl0 for x in x_t])
a_t = (Bl_vetor/m)*i_t - (k/m)*x_t - (b/m)*v_t

fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharex=True)
axs[0].plot(t_audio, Vin_amostras); axs[0].set_ylabel("Vin (V)")
axs[1].plot(t_audio, i_t);          axs[1].set_ylabel("i (A)")
axs[2].plot(t_audio, x_t);          axs[2].set_ylabel("x (m)")
axs[3].plot(t_audio, a_t);          axs[3].set_ylabel("ẍ (m/s²)")
axs[3].set_xlabel("tempo (s)")
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, f"resposta_temporal_{label_modelo}.png"))

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
plt.savefig(os.path.join(pasta_saida, f"espectros_{label_modelo}.png"))

a_norm = a_t / np.max(np.abs(a_t))
a_int16 = (a_norm * 32767).astype(np.int16)

nome_saida_audio = "TC02-out1.wav" if not NAO_LINEAR else "TC02-out2.wav"
wavfile.write(os.path.join(pasta_saida, nome_saida_audio), fs, a_int16)

fmin = 20
fmax = 22e3
npoints = 200
f_frf = np.logspace(np.log10(fmin), np.log10(fmax), npoints)
omega = 2 * np.pi * f_frf

A_lin = np.array([[-R/L, 0, -Bl0/L],
                   [0, 0, 1],
                   [Bl0/m, -k/m, -b/m]])
B_lin = np.array([1/L, 0, 0])
C_lin = np.array([0, 0, 1])
I3 = np.eye(3)

G = 1j * np.zeros(npoints)
for idx in range(npoints):
    aux = np.linalg.inv(1j*omega[idx]*I3 - A_lin)
    G[idx] = C_lin.dot(aux.dot(B_lin))

FRF = 20*np.log10(np.abs(1j*omega*G))
FRF = FRF - np.max(FRF)

Band_indexs = np.flatnonzero(np.where(FRF > -3, 1, 0))
fc_min = f_frf[Band_indexs[0]]
fc_max = f_frf[Band_indexs[-1]]
BW = fc_max - fc_min

print(f"Banda passante do alto-falante: {fc_min:.1f} Hz a {fc_max:.1f} Hz (largura: {BW:.1f} Hz)")

fig, ax1 = plt.subplots(figsize=(10, 5))

ax1.semilogx(f_vin, mag_vin / np.max(mag_vin), 'b', label="Espectro de Vin (normalizado)")
ax1.set_xlabel("Frequência (Hz)")
ax1.set_ylabel("Magnitude normalizada (áudio)", color='b')
ax1.set_xlim(10, 2e4)

ax2 = ax1.twinx()
ax2.semilogx(f_frf, FRF, 'g-D', markevery=[Band_indexs[0], Band_indexs[-1]], label="FRF do alto-falante")
ax2.axhline(-3, color='k', linestyle='-.', label="-3dB")
ax2.set_ylabel("FRF (dB)", color='g')

fig.legend(loc="upper right", bbox_to_anchor=(0.9, 0.9))
plt.title("Espectro do áudio vs. Banda passante do alto-falante")
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "frf_vs_espectro_audio.png"))
