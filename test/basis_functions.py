import matplotlib.pyplot as plt
import numpy as np

from finite_elem.equation import SturmLiouville


def p(x):
    return -1


def p_prime(x):
    return 0


def r(x):
    return 0


def f(x):
    return np.sin(np.pi * x)
    # return 12*x**2*np.sin(x) + 8*x**3*np.cos(x) - x**4*np.sin(x)
    # return -6*x**4+6*x**3+6*x**2-2


def ana_solution(x):
    return -np.sin(np.pi * x) / np.pi**2
    # return -x**3+x**2+x-1


a = -1
b = 1
N = 12
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
fm,grid = eq.fem()
xs = np.linspace(a, b, 200)
i = 1

us = np.array([eq.phi_quad(i, x) for x in xs])
us2 = np.array([eq.phi(i, x) for x in xs])
us3 = np.array([eq.phi_quart(i, x) for x in xs])

print(f'Linear = {eq.phi(i, grid[i])}; Quad {eq.phi_quad(i, grid[i])}, Quart {eq.phi_quart(i, grid[i])}')

plt.plot(grid,[0]*len(grid),"o",label="Grid points")
plt.plot(xs, us, label="Quad basis")
plt.plot(xs, us2, label="Linear basis")
plt.plot(xs, us3, label="Quartic basis")

plt.legend()
plt.grid()
plt.show()
