import os
import numpy as np
import matplotlib.pyplot as plt

pasta_script = os.path.dirname(os.path.abspath(__file__))
pasta_saida = os.path.join(pasta_script, "resultados")
os.makedirs(pasta_saida, exist_ok=True)

# 1. SOLUÇÃO MANUFATURADA

def u_exata(x):
    return 1 + 2*x - 3*x**2 + x**3 + 0.5*np.sin(2*np.pi*x) + np.exp(x/2)


def f(x):
    return -6 + 6*x - 2*np.pi**2*np.sin(2*np.pi*x) + 0.25*np.exp(x/2)

# 2. QUADRATURA DE GAUSS

def quadratura_gauss(a, b, ordem=5):
    xi, wi = np.polynomial.legendre.leggauss(ordem)

    x = (a + b)/2 + (b - a)/2 * xi
    w = (b - a)/2 * wi

    return x, w

# 3. MÉTODO DOS ELEMENTOS FINITOS - FEM

def fem(n_elem):
    n_nodes = n_elem + 1

    x = np.linspace(0.0, 1.0, n_nodes)

    K = np.zeros((n_nodes, n_nodes))
    F = np.zeros(n_nodes)

    # Montagem elemento por elemento
    for e in range(n_elem):

        x1 = x[e]
        x2 = x[e + 1]

        h = x2 - x1

        # Matriz de rigidez local
        Ke = (1.0 / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])

        # Quadratura para o vetor de força
        xg, wg = quadratura_gauss(x1, x2)

        N1 = (x2 - xg) / h
        N2 = (xg - x1) / h

        Fe = np.zeros(2)

        # Como u'' = f:
        #
        # integral(w' u') dx = - integral(w f) dx
        Fe[0] = -np.sum(wg * N1 * f(xg))
        Fe[1] = -np.sum(wg * N2 * f(xg))

        indices = [e, e + 1]

        # Montagem global
        for i_local in range(2):

            i_global = indices[i_local]

            F[i_global] += Fe[i_local]

            for j_local in range(2):

                j_global = indices[j_local]

                K[i_global, j_global] += Ke[i_local, j_local]

    # Condições de contorno de Dirichlet
    u = np.zeros(n_nodes)

    u[0] = u_exata(0.0)
    u[-1] = u_exata(1.0)

    # Sistema apenas para os nós internos
    K_interno = K[1:-1, 1:-1]

    F_interno = F[1:-1] - K[1:-1, 0]*u[0] - K[1:-1, -1]*u[-1]

    u[1:-1] = np.linalg.solve(K_interno, F_interno)

    return x, u


# ============================================================
# 4. MÉTODO DAS DIFERENÇAS FINITAS - FDM
# ============================================================

def fdm(n_elem):
    n_nodes = n_elem + 1

    x = np.linspace(0.0, 1.0, n_nodes)

    h = x[1] - x[0]

    u = np.zeros(n_nodes)

    u[0] = u_exata(0.0)
    u[-1] = u_exata(1.0)

    # Número de nós internos
    m = n_nodes - 2

    A = np.zeros((m, m))
    b = h**2 * f(x[1:-1])

    # u_{i-1} - 2u_i + u_{i+1} = h² f(x_i)
    for i in range(m):

        A[i, i] = -2.0

        if i > 0:
            A[i, i - 1] = 1.0

        if i < m - 1:
            A[i, i + 1] = 1.0

    # Inserção das condições de contorno
    b[0] -= u[0]
    b[-1] -= u[-1]

    u[1:-1] = np.linalg.solve(A, b)

    return x, u


# ============================================================
# 5. ERRO MÉDIO QUADRÁTICO
# ============================================================

def erro_medio_quadratico(x_nos, u_numerica, quantidade_pontos=5000):
    x = np.linspace(0.0, 1.0, quantidade_pontos)

    u_ref = u_exata(x)
    u_num = np.interp(x, x_nos, u_numerica)

    return np.mean((u_num - u_ref)**2)


# ============================================================
# 6. GRÁFICOS DA SOLUÇÃO MANUFATURADA E DO TERMO FONTE
# ============================================================

x_plot = np.linspace(0.0, 1.0, 1000)

plt.figure(figsize=(8, 5))

plt.plot(x_plot, u_exata(x_plot), label="Solução manufaturada")

plt.xlabel("x")
plt.ylabel("uₑ(x)")
plt.title("Solução manufaturada")

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "solucao_manufaturada.png"), dpi=300)
plt.show()
plt.figure(figsize=(8, 5))
plt.plot(x_plot, f(x_plot), label="f(x) = uₑ''(x)")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Termo fonte")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "termo_fonte.png"), dpi=300)
plt.show()

# 7. COMPARAÇÃO PARA n = 6 ELEMENTOS
n = 6
x_fem, u_fem = fem(n)
x_fdm, u_fdm = fdm(n)

x_exato = np.linspace(0.0, 1.0, 1000)
u_ref = u_exata(x_exato)


plt.figure(figsize=(8, 5))
plt.plot(x_exato, u_ref, label="Solução analítica")
plt.plot(x_fem, u_fem, "o--", label="FEM")
plt.plot(x_fdm, u_fdm, "s--", label="FDM")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Comparação entre solução analítica, FEM e FDM - n = 6")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "comparacao_n6.png"), dpi=300)
plt.show()

# 8. ANÁLISE DE CONVERGÊNCIA

valores_n = [4, 8, 16, 32]

erros_fem = []
erros_fdm = []

for n in valores_n:

    x_fem, u_fem = fem(n)
    x_fdm, u_fdm = fdm(n)

    erros_fem.append(erro_medio_quadratico(x_fem, u_fem))
    erros_fdm.append(erro_medio_quadratico(x_fdm, u_fdm))

# 9. ORDEM DE CONVERGÊNCIA

rmse_fem = np.sqrt(erros_fem)
rmse_fdm = np.sqrt(erros_fdm)

ordem_fem = [np.nan]
ordem_fdm = [np.nan]

for i in range(1, len(valores_n)):

    p_fem = np.log(rmse_fem[i - 1] / rmse_fem[i]) / np.log(2)
    p_fdm = np.log(rmse_fdm[i - 1] / rmse_fdm[i]) / np.log(2)

    ordem_fem.append(p_fem)
    ordem_fdm.append(p_fdm)

# 10. TABELA DE RESULTADOS

print() 
print("analise de convergencia:")
print(f"{'n':>5} {'EMQ FEM':>16} {'Ordem FEM':>12} {'EMQ FDM':>16} {'Ordem FDM':>12}")

for i, n in enumerate(valores_n):

    if i == 0:

        print(
            f"{n:5d} "
            f"{erros_fem[i]:16.8e} "
            f"{'-':>12} "
            f"{erros_fdm[i]:16.8e} "
            f"{'-':>12}"
        )

    else:

        print(
            f"{n:5d} "
            f"{erros_fem[i]:16.8e} "
            f"{ordem_fem[i]:12.4f} "
            f"{erros_fdm[i]:16.8e} "
            f"{ordem_fdm[i]:12.4f}"
        )

# 11. GRÁFICO DE CONVERGÊNCIA

plt.figure(figsize=(8, 5))
plt.loglog(valores_n, erros_fem, "o-", label="FEM")
plt.loglog(valores_n, erros_fdm, "s-", label="FDM")
plt.xlabel("Número de elementos (n)")
plt.ylabel("Erro médio quadrático")
plt.title("Convergência dos métodos FEM e FDM")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "convergencia_fem_fdm.png"), dpi=300)
plt.show()