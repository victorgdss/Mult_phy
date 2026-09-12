import os
import numpy as np
import matplotlib.pyplot as plt

pasta_script = os.path.dirname(os.path.abspath(__file__))
pasta_saida = os.path.join(pasta_script, "resultados")
os.makedirs(pasta_saida, exist_ok=True)
# 1. DADOS DO PROBLEMA
L = 1.0
Va = 10.0
Vb = 0.0
eps1 = 3.0
eps2 = 1.0

# 2. PERMISSIVIDADE RELATIVA

def epsilon_r(z):
    z = np.asarray(z)

    return np.where(z < L/2, eps1, eps2)


# 3. SOLUÇÃO ANALÍTICA

def phi_exata(z):
    z = np.asarray(z)

    return np.where(z <= L/2, 10.0 - 5.0*z, 15.0*(1.0 - z))


# 4. MÉTODO DOS ELEMENTOS FINITOS - FEM

def fem(n_nodes):
    if n_nodes % 2 == 0:
        raise ValueError("Na Atividade 2, n deve ser ímpar.")

    z = np.linspace(0.0, L, n_nodes)

    n_elem = n_nodes - 1

    K = np.zeros((n_nodes, n_nodes))
    F = np.zeros(n_nodes)

    # Montagem elemento por elemento
    for e in range(n_elem):

        z1 = z[e]
        z2 = z[e + 1]

        h = z2 - z1

        # Avaliação da permissividade no centro do elemento
        z_medio = (z1 + z2) / 2

        eps = float(epsilon_r(z_medio))

        # Matriz local
        Ke = (eps / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])

        indices = [e, e + 1]

        # Montagem global
        for i_local in range(2):

            i_global = indices[i_local]

            for j_local in range(2):

                j_global = indices[j_local]

                K[i_global, j_global] += Ke[i_local, j_local]

    # Condições de contorno
    phi = np.zeros(n_nodes)

    phi[0] = Va
    phi[-1] = Vb

    # Sistema interno
    K_interno = K[1:-1, 1:-1]

    F_interno = F[1:-1] - K[1:-1, 0]*phi[0] - K[1:-1, -1]*phi[-1]

    phi[1:-1] = np.linalg.solve(K_interno, F_interno)

    return z, phi


# 5. MÉTODO DAS DIFERENÇAS FINITAS - FDM CONSERVATIVO

def fdm(n_nodes):

    if n_nodes % 2 == 0:
        raise ValueError("Na Atividade 2, n deve ser ímpar.")

    z = np.linspace(0.0, L, n_nodes)

    n_internos = n_nodes - 2

    A = np.zeros((n_internos, n_internos))
    b = np.zeros(n_internos)

    for i in range(1, n_nodes - 1):

        linha = i - 1

        # Permissividade no intervalo à esquerda
        z_esq = (z[i - 1] + z[i]) / 2

        # Permissividade no intervalo à direita
        z_dir = (z[i] + z[i + 1]) / 2

        eps_esq = float(epsilon_r(z_esq))
        eps_dir = float(epsilon_r(z_dir))

        # Termo central
        A[linha, linha] = -(eps_esq + eps_dir)

        # Nó à esquerda
        if linha > 0:
            A[linha, linha - 1] = eps_esq
        else:
            b[linha] -= eps_esq * Va

        # Nó à direita
        if linha < n_internos - 1:
            A[linha, linha + 1] = eps_dir
        else:
            b[linha] -= eps_dir * Vb

    phi = np.zeros(n_nodes)

    phi[0] = Va
    phi[-1] = Vb

    phi[1:-1] = np.linalg.solve(A, b)

    return z, phi

# 6. SOLUÇÃO NUMÉRICA

n = 15

z_fem, phi_fem = fem(n)
z_fdm, phi_fdm = fdm(n)

phi_ref_nos = phi_exata(z_fem)

# 7. ERROS NOS NÓS

erro_fem = np.max(np.abs(phi_fem - phi_ref_nos))
erro_fdm = np.max(np.abs(phi_fdm - phi_ref_nos))

print()
print("Capacitor com dois dieletricos")
print(f"Número de nós: {n}")
print(f"Número de elementos: {n - 1}")
print()

print(f"Erro máximo FEM = {erro_fem:.10e}")
print(f"Erro máximo FDM = {erro_fdm:.10e}")

# 8. TABELA DOS RESULTADOS

print()
print(f"{'z':>8} {'Exata':>14} {'FEM':>14} {'FDM':>14}")

print("-" * 54)

for i in range(n):

    print(
        f"{z_fem[i]:8.4f} "
        f"{phi_ref_nos[i]:14.8f} "
        f"{phi_fem[i]:14.8f} "
        f"{phi_fdm[i]:14.8f}"
    )

# 9. INTERFACE DOS DIELÉTRICOS
indice_interface = np.argmin(np.abs(z_fem - L/2))

print()

print("INTERFACE:")
print(f"Posição da interface = {z_fem[indice_interface]:.6f}")
print(f"Solução exata        = {phi_ref_nos[indice_interface]:.8f} V")
print(f"FEM                  = {phi_fem[indice_interface]:.8f} V")
print(f"FDM                  = {phi_fdm[indice_interface]:.8f} V")



# 10. VERIFICAÇÃO DO FLUXO NA INTERFACE

i = indice_interface

h_esq = z_fem[i] - z_fem[i - 1]
h_dir = z_fem[i + 1] - z_fem[i]

dphi_esq = (phi_fem[i] - phi_fem[i - 1]) / h_esq
dphi_dir = (phi_fem[i + 1] - phi_fem[i]) / h_dir

fluxo_esq = eps1 * dphi_esq
fluxo_dir = eps2 * dphi_dir

print()
print("continuidade do fluxo - FEM:")
print(f"epsilon1 * phi'(L/2-) = {fluxo_esq:.8f}")
print(f"epsilon2 * phi'(L/2+) = {fluxo_dir:.8f}")

print(f"Diferença = {abs(fluxo_esq - fluxo_dir):.10e}")


# 11. DIFERENÇA FEM x FDM

diferenca_fem_fdm = np.max(np.abs(phi_fem - phi_fdm))

print()
print(f"Máxima diferença FEM x FDM = {diferenca_fem_fdm:.10e}")

# 12. GRÁFICO COMPARATIVO


z_exato = np.linspace(0.0, L, 1000)

plt.figure(figsize=(8, 5))
plt.plot(z_exato, phi_exata(z_exato), label="Solução analítica")
plt.plot(z_fem, phi_fem, "o-", label="FEM")
plt.plot(z_fdm, phi_fdm, "s--", label="FDM conservativo")
plt.axvline(L/2, linestyle=":", label="Interface z = L/2")
plt.xlabel("z")
plt.ylabel("Potencial elétrico φ(z) [V]")
plt.title("Capacitor composto por dois dielétricos")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "comparacao_atividade2.png"))
plt.show()