import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Parâmetros dE entrada
m    = 0.073
g    = 9.81
k    = 6.51e-5
R    = 11.0
imax = 2.0
x0   = 8.5e-3

#Funções das equações físicas
def fm(i, x):
    return k/2 * (i**2 / x**2)

def Bl(i, x):
    return k/2 * (i / x**2)

def L(x):
    return k / x

i0 = x0 * np.sqrt(2 * m * g / k)
u0 = R * i0
print("Ponto de equilibrio:")
print(f"  i0 = {i0:.4f} A")
print(f"  u0 = {u0:.4f} V")

#Parâmetros do controlador PID
Kp = 4000
Ki = 800
Kd = 100

#Funcao do sistema com PID
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

def simula_epilson(eps):
    y_0 = [i0, x0 + eps, 0.0, 0.0]
    sol = solve_ivp(f_pid, [0.0, 2.0], y_0, max_step=1e-3)
    i_pico = np.max(sol.y[0,:])
    desvio_final = abs(sol.y[1,-1] - x0)
    estavel = (i_pico < imax) and (desvio_final < 1e-5)

    return estavel, i_pico

# ------------------------------------------------------------
# Demonstre que o sistema é estável para um dado valor de ϵ
# ------------------------------------------------------------
epilson_teste = 0.0002  
y_0 = [i0, x0 + epilson_teste, 0.0, 0.0]
sol = solve_ivp(f_pid, [0.0, 2.0], y_0, max_step=1e-3)
fig = plt.figure(figsize=(15, 4))
fig.suptitle(f'Parte 1 - Estabilidade para epilson = {epilson_teste*1000:.1f} mm')
plt.subplot(1, 3, 1)
plt.plot(sol.t, sol.y[1,:] - x0, 'k', linewidth=2)
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('x(t) - x0  (m)')
plt.title('Posicao (desvio)')

plt.subplot(1, 3, 2)
plt.plot(sol.t, sol.y[2,:], 'k', linewidth=2)
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('v(t)  (m/s)')
plt.title('Velocidade')

plt.subplot(1, 3, 3)
plt.plot(sol.t, sol.y[0,:], 'k', linewidth=2)
plt.axhline(imax, color='r', linestyle='--', label='i_max = 2A')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('i(t)  (A)')
plt.title('Corrente')
plt.legend()

plt.tight_layout()
plt.savefig('parte1_demonstracao.png')
print("Grafico salvo como parte1_demonstracao.png")

# ----------------------------------------------------------------------------------------------------------------------------
# Avalie qual a faixa de valores de ϵ na qual o sistema se mantem estável. Considere a corrente máxima da fonte igual a 2A.
# ----------------------------------------------------------------------------------------------------------------------------
print("\nVarredura grossa:")
epilsons_grosso = np.linspace(-0.006, 0.004, 50)
estavel_grosso = []
i_pico_grosso = []
encontrou_estavel = False
idx_primeiro      = None
idx_ultimo        = None

for j, eps in enumerate(epilsons_grosso):

    estavel, i_pico = simula_epilson(eps)
    estavel_grosso.append(estavel)
    i_pico_grosso.append(i_pico)
    print(f"  eps = {eps*1000:+.3f} mm  /  i_pico = {i_pico:.4f} A  |  {'ESTAVEL' if estavel else 'INSTAVEL'}")
    if estavel:
        encontrou_estavel = True
        if idx_primeiro is None:
            idx_primeiro = j
        idx_ultimo = j

    if encontrou_estavel and not estavel:
        break

print(f"\n  Transicao inferior: entre {epilsons_grosso[idx_primeiro-1]*1000:+.3f} mm e {epilsons_grosso[idx_primeiro]*1000:+.3f} mm")
print(f"  Transicao superior: entre {epilsons_grosso[idx_ultimo]*1000:+.3f} mm e {epilsons_grosso[idx_ultimo+1]*1000:+.3f} mm")

print("\nVarredura fina:")
epilsons_finos_inf = np.linspace(epilsons_grosso[idx_primeiro-1], epilsons_grosso[idx_primeiro], 100)
limite_inferior = epilsons_grosso[idx_primeiro]
for eps in epilsons_finos_inf:
    estavel, _ = simula_epilson(eps)
    if estavel:
        limite_inferior = eps
        break
epilsons_finos_sup = np.linspace(epilsons_grosso[idx_ultimo], epilsons_grosso[idx_ultimo+1], 100)
limite_superior = epilsons_grosso[idx_ultimo]
for eps in epilsons_finos_sup:
    estavel, _ = simula_epilson(eps)
    if not estavel:
        break
    limite_superior = eps
print(f"\n Limite inferior: epilson = {limite_inferior*1000:+.4f} mm")
print(f"  Limite superior: epilson = {limite_superior*1000:+.4f} mm")
print(f"  Faixa estavel:   {(limite_superior - limite_inferior)*1000:.4f} mm")


n = len(i_pico_grosso)
fig = plt.figure(figsize=(10, 5))
fig.suptitle(f'Parte 2 - Corrente de pico vs epilson  |  Kp={Kp}, Ki={Ki}, Kd={Kd}')
plt.subplot(1, 1, 1)
plt.plot(epilsons_grosso[:n]*1000, i_pico_grosso, 'k', linewidth=2, label='Corrente de pico')
plt.axhline(imax, color='r', linestyle='--', label='i_max = 2A')
plt.axvline(limite_inferior*1000, color='b', linestyle='--', label=f'eps_inf = {limite_inferior*1000:.2f} mm')
plt.axvline(limite_superior*1000, color='g', linestyle='--', label=f'eps_sup = {limite_superior*1000:.2f} mm')
plt.grid(True)
plt.xlabel('epilson (mm)')
plt.ylabel('i_pico (A)')
plt.title('Corrente de pico em funcao de epilson')
plt.legend()

plt.tight_layout()
plt.savefig('parte2_faixa_epilson.png')