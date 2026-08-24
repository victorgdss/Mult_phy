import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
np.random.seed(42)

# Parametros do TP
m  = 0.073
g = 9.81
k = 6.51e-5
R = 11.0
imax = 2.0
x0 = 8.5e-3

# Funcoes
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

#Ganhos do controlador
Kp = 4000
Ki = 800
Kd = 100
dp = 0.04 * x0


#Sistema com ruido na medicao de posicao
def f_pid_ruido(t, y):
    i = y[0]
    x = y[1]
    v = y[2]
    E = y[3]
    ruido = np.random.normal(0, dp)
    x_med = x + ruido
    Lx = L(x)
    fmix = fm(i, x)
    Blix = Bl(i, x)
    e_med = x_med - x0
    u = u0 + Kp*e_med + Ki*E + Kd*v
    di = u/Lx - R/Lx*i - Blix/Lx*v
    dx = v
    dv = g - fmix/m
    dE = e_med

    return [di, dx, dv, dE]


# ------------------------------------------------------------
# Reavalie a faixa de valores de ϵ na qual o sistema é estável.
# ------------------------------------------------------------
epilson_teste = 0.0002    
y_0 = [i0, x0 + epilson_teste, 0.0, 0.0]
sol = solve_ivp(f_pid_ruido, [0.0, 2.0], y_0, max_step=1e-3)

fig = plt.figure(figsize=(15, 4))
fig.suptitle('Parte 1 - Estabilidade com ruido (epilson = 0.2 mm, dp = 0.3400 mm)')
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
plt.savefig('atividade2_parte1.png')


def simula_epilson(eps):
    y_0 = [i0, x0 + eps, 0.0, 0.0]
    sol = solve_ivp(f_pid_ruido, [0.0, 5.0], y_0, max_step=1e-3)
    i_pico = np.max(sol.y[0,:])
    desvio_final = abs(sol.y[1,-1] - x0)
    estavel = (i_pico < imax) and (desvio_final < 1e-3)

    return estavel, i_pico

print("\nVarredura grossa de epilson :")

epilsons_grosso = np.linspace(-0.006, 0.004, 50)
estavel_grosso = []
i_pico_grosso = []
idx_primeiro_estavel = None
idx_ultimo_estavel   = None
viu_estavel = False

for idx, eps in enumerate(epilsons_grosso):
    estavel, i_pico = simula_epilson(eps)
    estavel_grosso.append(estavel)
    i_pico_grosso.append(i_pico)
    print(f"  epilson = {eps*1000:+.3f} mm  /  i_pico = {i_pico:.4f} A  |  {'ESTAVEL' if estavel else 'INSTAVEL'}")

    if estavel and idx_primeiro_estavel is None:
        idx_primeiro_estavel = idx

    if estavel:
        viu_estavel = True
        idx_ultimo_estavel = idx
    elif viu_estavel:
        break    

print(f"\n  Transicao inferior: entre {epilsons_grosso[idx_primeiro_estavel-1]*1000:+.3f} mm e {epilsons_grosso[idx_primeiro_estavel]*1000:+.3f} mm")
print(f"  Transicao superior: entre {epilsons_grosso[idx_ultimo_estavel]*1000:+.3f} mm e {epilsons_grosso[idx_ultimo_estavel+1]*1000:+.3f} mm")


print("\nVarredura fina nos limites:")

epilsons_finos_inf = np.linspace(epilsons_grosso[idx_primeiro_estavel - 1], epilsons_grosso[idx_primeiro_estavel], 100)
limite_inferior = epilsons_grosso[idx_primeiro_estavel]
for eps in epilsons_finos_inf:
    estavel, _ = simula_epilson(eps)
    if estavel:
        limite_inferior = eps
        break

epilsons_finos_sup = np.linspace(epilsons_grosso[idx_ultimo_estavel], epilsons_grosso[idx_ultimo_estavel + 1], 100)
limite_superior = epilsons_grosso[idx_ultimo_estavel]
for eps in epilsons_finos_sup[::-1]:
    estavel, _ = simula_epilson(eps)
    if estavel:
        limite_superior = eps
        break

print(f"\n  Limite inferior: epilson = {limite_inferior*1000:+.4f} mm")
print(f"  Limite superior: epilson = {limite_superior*1000:+.4f} mm")
print(f"  Faixa estavel total: {(limite_superior - limite_inferior)*1000:.4f} mm")


epilsons_mm = epilsons_grosso[:len(i_pico_grosso)] * 1000
i_pico_plot = np.clip(np.array(i_pico_grosso), 0, 2.5)
fig = plt.figure(figsize=(10, 5))
fig.suptitle('Parte 2 - Faixa de epilson estavel com ruido')
plt.subplot(1, 1, 1)
plt.axvspan(limite_inferior*1000, limite_superior*1000, alpha=0.15, color='green', label='Regiao estavel')
plt.plot(epilsons_mm, i_pico_plot, 'k', linewidth=2, label='Corrente de pico')
plt.axhline(imax, color='r', linestyle='--', linewidth=1.5, label='i_max = 2 A')
plt.axvline(limite_inferior*1000, color='b', linestyle='--', linewidth=1.5, label=f'epilson_inf = {limite_inferior*1000:.2f} mm')
plt.axvline(limite_superior*1000, color='g', linestyle='--', linewidth=1.5, label=f'epilson_sup = {limite_superior*1000:.2f} mm')
plt.xlabel('epilson (mm)')
plt.ylabel('Corrente de pico (A)')
plt.title('Corrente de pico vs epilson')
plt.ylim([0, 2.5])
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('atividade2_parte2.png')

# -----------------------------------------------------------------------
# Para um ϵ constante, avalie a estabilidade do sistema em relacaoo a sd.
# -----------------------------------------------------------------------
epilson_parte3 = 0.0018
print("\nVarredura de sd com epilson fixo em 1.8 mm:")
sds = np.linspace(0, 20 * dp, 30)
estavel_sd = []
i_pico_sd  = []

for sd_teste in sds:

    def funcao_com_ruido(t, y):
        i = y[0]
        x = y[1]
        v = y[2]
        E = y[3]

        ruido = np.random.normal(0, sd_teste)
        x_med = x + ruido
        Lx = L(x)
        fmix = fm(i, x)
        Blix = Bl(i, x)
        e_med = x_med - x0
        u = u0 + Kp*e_med + Ki*E + Kd*v
        di = u/Lx - R/Lx*i - Blix/Lx*v
        dx = v
        dv = g - fmix/m
        dE = e_med

        return [di, dx, dv, dE]

    y_0 = [i0, x0 + epilson_parte3, 0.0, 0.0]
    sol = solve_ivp(funcao_com_ruido, [0.0, 1.0], y_0, max_step=5e-3)
    i_pico = np.max(sol.y[0,:])
    desvio_final = abs(sol.y[1,-1] - x0)
    estavel = (i_pico < imax) and (desvio_final < 1e-3)
    estavel_sd.append(estavel)
    i_pico_sd.append(i_pico)

    print(f"  sd = {sd_teste*1000:.4f} mm  ({sd_teste/dp:.2f}x dp)  /  i_pico = {i_pico:.4f} A  |  {'ESTAVEL' if estavel else 'INSTAVEL'}")

    if not estavel:
        break    

ultimo_sd_estavel = 0
for i, est in enumerate(estavel_sd):
    if est:
        ultimo_sd_estavel = sds[i]

sds_mm = sds[:len(i_pico_sd)] * 1000
i_pico_sd_plot = np.clip(np.array(i_pico_sd), 0, 2.5)

fig = plt.figure(figsize=(10, 5))
fig.suptitle('Parte 3 - Estabilidade vs desvio padrao do ruido (epilson = 1.8 mm)')

plt.subplot(1, 1, 1)
plt.axvspan(0, ultimo_sd_estavel*1000,alpha=0.15, color='green', label='Regiao estavel')
plt.plot(sds_mm, i_pico_sd_plot, 'k', linewidth=2, label='Corrente de pico')
plt.axhline(imax, color='r', linestyle='--', linewidth=1.5, label='i_max = 2 A')
plt.axvline(dp*1000, color='b', linestyle=':', linewidth=1.5, label=f'sd nominal = {dp*1000:.4f} mm')
plt.xlabel('sd (mm)')
plt.ylabel('Corrente de pico (A)')
plt.title('Corrente de pico vs sd')
plt.ylim([0, 2.5])
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('atividade2_parte3.png')