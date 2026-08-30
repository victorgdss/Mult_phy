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

freq_sinal = 150
amplitude = 2.0
n_periodos = 30
fs = 44100
duracao = n_periodos / freq_sinal

t_audio = np.arange(0, duracao, 1/fs)
Vin_amostras = amplitude * np.sin(2*np.pi*freq_sinal*t_audio)

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

def derivadas(t, z, nao_linear):
    i, x, v = z
    u = u_interp(t)
    Bl = Bl_de_x(x) if nao_linear else Bl0

    di_dt = -(R/L)*i - (Bl/L)*v + u/L
    dx_dt = v
    dv_dt = (Bl/m)*i - (k/m)*x - (b/m)*v

    return [di_dt, dx_dt, dv_dt]

z0 = [0.0, 0.0, 0.0]
t_span = (t_audio[0], t_audio[-1])

sol_lin = solve_ivp(derivadas, t_span, z0, args=(False,),
                     t_eval=t_audio, method="RK45", max_step=1/fs)

sol_nlin = solve_ivp(derivadas, t_span, z0, args=(True,),
                      t_eval=t_audio, method="RK45", max_step=1/fs)

def calc_aceleracao(sol, nao_linear):
    i_t, x_t, v_t = sol.y[0], sol.y[1], sol.y[2]
    Bl_vetor = np.array(
        [Bl_de_x(x) if nao_linear else Bl0 for x in x_t]
    )
    a_t = (Bl_vetor/m)*i_t - (k/m)*x_t - (b/m)*v_t
    return i_t, x_t, v_t, a_t

i_lin, x_lin, v_lin, a_lin = calc_aceleracao(sol_lin, False)
i_nlin, x_nlin, v_nlin, a_nlin = calc_aceleracao(sol_nlin, True)

print(f"Pico de x (linear):     {np.max(np.abs(x_lin)):.6e} m")
print(f"Pico de x (não linear): {np.max(np.abs(x_nlin)):.6e} m")

fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

axs[0].plot(t_audio, Vin_amostras, 'k')
axs[0].set_ylabel("Vin (V)")
axs[0].set_title(f"Sinal senoidal de entrada ({freq_sinal} Hz, {amplitude}V pico)")

axs[1].plot(t_audio, x_lin, label="Linear")
axs[1].plot(t_audio, x_nlin, label="Não linear", alpha=0.7)
axs[1].axhline(x1, color='gray', linestyle='--', linewidth=0.8)
axs[1].axhline(-x1, color='gray', linestyle='--', linewidth=0.8)
axs[1].set_ylabel("x (m)")
axs[1].legend()
axs[1].set_title("Posição do cone: linear vs não linear (linhas tracejadas = x1)")

axs[2].plot(t_audio, a_lin, label="Linear")
axs[2].plot(t_audio, a_nlin, label="Não linear", alpha=0.7)
axs[2].set_ylabel("ẍ (m/s²)")
axs[2].set_xlabel("tempo (s)")
axs[2].legend()
axs[2].set_title("Aceleração do cone: linear vs não linear")

plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "atividade4_comparacao_senoide.png"))

def espectro(sinal, fs):
    N = len(sinal)
    freqs = np.fft.rfftfreq(N, d=1/fs)
    mag = np.abs(np.fft.rfft(sinal)) / N
    return freqs, mag

f_vin, mag_vin = espectro(Vin_amostras, fs)
f_lin, mag_lin = espectro(a_lin, fs)
f_nlin, mag_nlin = espectro(a_nlin, fs)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(f_vin, mag_vin / np.max(mag_vin), label="Vin (entrada)")
ax.plot(f_lin, mag_lin / np.max(mag_lin), label="ẍ linear")
ax.plot(f_nlin, mag_nlin / np.max(mag_nlin), '--', label="ẍ não linear")
ax.set_xlim(0, 8 * freq_sinal)
ax.set_xlabel("Frequência (Hz)")
ax.set_ylabel("Magnitude normalizada")
ax.set_title("Espectro de frequência: entrada vs. saída (linear e não linear)")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "atividade4_espectros.png"))

n_periodos_zoom = 3
t_zoom_final = duracao
t_zoom_inicio = t_zoom_final - n_periodos_zoom/freq_sinal
mask = (t_audio >= t_zoom_inicio) & (t_audio <= t_zoom_final)

fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
axs[0].plot(t_audio[mask], x_lin[mask], label="Linear")
axs[0].plot(t_audio[mask], x_nlin[mask], label="Não linear")
axs[0].set_ylabel("x (m)")
axs[0].legend()
axs[0].set_title("Regime permanente (zoom) - Posição")

axs[1].plot(t_audio[mask], a_lin[mask], label="Linear")
axs[1].plot(t_audio[mask], a_nlin[mask], label="Não linear")
axs[1].set_ylabel("ẍ (m/s²)")
axs[1].set_xlabel("tempo (s)")
axs[1].legend()
axs[1].set_title("Regime permanente (zoom) - Aceleração")

plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "atividade4_zoom_regime_permanente.png"))

n_periodos_zoom_ini = 3
t_zoom_ini_fim = n_periodos_zoom_ini / freq_sinal
mask_ini = t_audio <= t_zoom_ini_fim

fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
axs[0].plot(t_audio[mask_ini], x_lin[mask_ini], label="Linear")
axs[0].plot(t_audio[mask_ini], x_nlin[mask_ini], label="Não linear")
axs[0].set_ylabel("x (m)")
axs[0].legend()
axs[0].set_title("Regime transitório (zoom) - Posição")

axs[1].plot(t_audio[mask_ini], a_lin[mask_ini], label="Linear")
axs[1].plot(t_audio[mask_ini], a_nlin[mask_ini], label="Não linear")
axs[1].set_ylabel("ẍ (m/s²)")
axs[1].set_xlabel("tempo (s)")
axs[1].legend()
axs[1].set_title("Regime transitório (zoom) - Aceleração")

plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "atividade4_zoom_transitorio.png"))
