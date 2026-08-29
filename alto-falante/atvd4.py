import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# =========================================================
# 1) PARÂMETROS DO ALTO-FALANTE (iguais aos anteriores)
# =========================================================

m  = 14.35e-3   # massa do cone+bobina (kg)
b  = 0.786      # amortecimento mecânico (kg/s)
k  = 1852       # rigidez da suspensão (N/m)
Bl0 = 4.95      # fator de força (N/A) em repouso, usado no modelo linear
L  = 266e-6     # indutância da bobina (H)
R  = 3.3        # resistência da bobina (ohm)

pasta_script = os.path.dirname(os.path.abspath(__file__))
pasta_saida = os.path.join(pasta_script, "resultados")
os.makedirs(pasta_saida, exist_ok=True)

# =========================================================
# 2) GERANDO O SINAL SENOIDAL DE ENTRADA
# =========================================================

freq_sinal = 150  # Hz - dentro da faixa de voz humana (85-255 Hz)
amplitude = 2.0   # V (pico) -> 4 Vpp, conforme pedido
n_periodos = 30   # mínimo 20 períodos exigido, uso 30 pra folga
fs = 44100        # taxa de amostragem (Hz) - padrão de áudio
duracao = n_periodos / freq_sinal   # duração total do sinal (s)

t_audio = np.arange(0, duracao, 1/fs)
Vin_amostras = amplitude * np.sin(2*np.pi*freq_sinal*t_audio)

u_interp = interp1d(t_audio, Vin_amostras, bounds_error=False, fill_value=0.0)

# =========================================================
# 3) MODELO NÃO LINEAR DE Bl(x) - mesmo da Atividade 2/3
# =========================================================

# Usamos o mesmo x_max obtido na Atividade 1 (simulação com áudio real).
# Se quiser, pode recalcular rodando este script primeiro em modo linear
# e pegando o x_max específico para a senoide - mas para manter a
# CONSISTÊNCIA do modelo não linear entre atividades, o ideal é usar o
# mesmo x1/x2 definidos anteriormente.

x_max_linear = 4.3514e-4   # valor obtido na Atividade 1 (modelo linear com áudio real)
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

# =========================================================
# 4) FUNÇÃO DE DERIVADAS (com parâmetro para escolher o modelo)
# =========================================================

def derivadas(t, z, nao_linear):
    i, x, v = z
    u = u_interp(t)
    Bl = Bl_de_x(x) if nao_linear else Bl0

    di_dt = -(R/L)*i - (Bl/L)*v + u/L
    dx_dt = v
    dv_dt = (Bl/m)*i - (k/m)*x - (b/m)*v

    return [di_dt, dx_dt, dv_dt]

# =========================================================
# 5) RESOLVER A ODE PARA OS DOIS MODELOS (linear e não linear)
# =========================================================

z0 = [0.0, 0.0, 0.0]
t_span = (t_audio[0], t_audio[-1])

sol_lin = solve_ivp(derivadas, t_span, z0, args=(False,),
                     t_eval=t_audio, method="RK45", max_step=1/fs)

sol_nlin = solve_ivp(derivadas, t_span, z0, args=(True,),
                      t_eval=t_audio, method="RK45", max_step=1/fs)

# Aceleração para cada modelo
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

# =========================================================
# 6) COMPARAÇÃO GRÁFICA: LINEAR x NÃO LINEAR
# =========================================================

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

# =========================================================
# 7) ZOOM NO REGIME PERMANENTE (últimos períodos, sem transitório)
# =========================================================

# O regime transitório ocorre nos primeiros períodos, até o sistema
# "assentar". Aqui damos zoom nos últimos períodos para ver o regime
# permanente já estabilizado.

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

print("Atividade 4 concluída.")
