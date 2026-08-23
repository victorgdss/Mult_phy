import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- Parâmetros físicos ---
m    = 0.073
g    = 9.81
k    = 6.51e-5
R    = 11.0
imax = 2.0
x0   = 8.5e-3

# --- Funções físicas ---
def fm(i, x):
    return k/2 * (i**2 / x**2)

def Bl(i, x):
    return k/2 * (i / x**2)

def L(x):
    return k / x

# --- Ponto de equilíbrio ---
i0 = x0 * np.sqrt(2 * m * g / k)
u0 = R * i0

# --- Conjuntos de ganhos a testar ---
# (Kp, Ki, Kd)
conjuntos = [
    (2000,  400,  50),
    (3000,  600,  80),
    (4000,  800, 100),
    (5000, 1000, 120),
    (6000, 1200, 150),
    (7000, 1500, 180),
    (8000, 2000, 200),
]

# --- Condição inicial ---
epsilon = 0.001
y_0 = [i0, x0 + epsilon, 0.0, 0.0]

t_0   = 0.0
t_end = 1.0

# --- Gráficos ---
fig = plt.figure(figsize=(15, 4))
fig.suptitle(f'Comparação de ganhos PID')

for Kp, Ki, Kd in conjuntos:

    def f_pid(t, y):
        i = y[0]
        x = y[1]
        v = y[2]
        E = y[3]

        Lx   = L(x)
        fmix = fm(i, x)
        Blix = Bl(i, x)

        e = x - x0
        u = u0 + Kp*e + Ki*E + Kd*v

        di = u/Lx - R/Lx*i - Blix/Lx*v
        dx = v
        dv = g - fmix/m
        dE = e

        return [di, dx, dv, dE]

    sol = solve_ivp(f_pid, [t_0, t_end], y_0, max_step=5e-3)

    label = f'Kp={Kp} Ki={Ki} Kd={Kd}'

    plt.subplot(1, 3, 1)
    plt.plot(sol.t, sol.y[1,:] - x0, linewidth=2, label=label)

    plt.subplot(1, 3, 2)
    plt.plot(sol.t, sol.y[2,:], linewidth=2, label=label)

    plt.subplot(1, 3, 3)
    plt.plot(sol.t, sol.y[0,:], linewidth=2, label=label)

# --- Finaliza os gráficos ---
plt.subplot(1, 3, 1)
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('x(t) - x0  (m)')
plt.title('Posição (desvio)')
plt.legend(fontsize=7)

plt.subplot(1, 3, 2)
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('v(t)  (m/s)')
plt.title('Velocidade')
plt.legend(fontsize=7)

plt.subplot(1, 3, 3)
plt.axhline(imax, color='r', linestyle='--', label='i_max = 2A')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('i(t)  (A)')
plt.title('Corrente')
plt.legend(fontsize=7)

plt.tight_layout()
plt.savefig('comparacao_pid.png')
print("Gráfico salvo como comparacao_pid.png")