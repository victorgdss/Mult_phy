import os
import numpy as np
import matplotlib.pyplot as plt

Bl0 = 4.95
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

x_vals = np.linspace(-1.3*x2, 1.3*x2, 1000)
Bl_vals = np.array([Bl_de_x(x) for x in x_vals])

pasta_script = os.path.dirname(os.path.abspath(__file__))
pasta_saida = os.path.join(pasta_script, "resultados")
os.makedirs(pasta_saida, exist_ok=True)

plt.figure(figsize=(8, 5))
plt.plot(x_vals, Bl_vals, 'b', linewidth=2)

for xv in [-x2, -x1, x1, x2]:
    plt.axvline(xv, color='gray', linestyle='--', linewidth=0.8)

plt.xlabel("x (m)")
plt.ylabel("Bl (N/A)")
plt.title("Fator de força não linear Bl(x)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(pasta_saida, "bl_de_x.png"))
plt.show()

print(f"x1 = {x1:.4e} m")
print(f"x2 = {x2:.4e} m")
