import matplotlib.pyplot as plt
import numpy as np

from finite_elem.equation import SturmLiouville


def p(x):
    return -11


def p_prime(x):
    return 0


def r(x):
    return 3 * x


def f(x):
    return 3 * np.sin(np.pi * x) * x**3 - 11 * (
        (2 - np.pi**2 * x**2) * np.sin(np.pi * x) + 4 * x * np.pi * np.cos(np.pi * x)
    )


def ana_solution(x):
    return -np.sin(np.pi * x) * x**2


a = -1
b = 1
N = 25
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
us, xs1 = eq.fem()

us2, xs = eq.sem()


plt.plot(
    np.linspace(a, b, N + 1),
    [ana_solution(x) for x in np.linspace(a, b, N + 1)],
    label="Analytical",
    alpha=0.5,
)
plt.plot(np.linspace(a, b, N + 1), us, label="FEM", linestyle=(0, (3, 10)))
plt.plot(xs, us2, label="SEM", linestyle=(1, (5, 10)))

plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
plt.grid()
plt.show()
plt.clf()


def p2(x):
    return 1


def p_prime2(x):
    return 0


def r2(x):
    return 6 * x


def f2(x):
    return -6 * x**4 + 6 * x**3 + 6 * x**2 - 2


def ana_solution2(x):
    return -(x**3) + x**2 + x - 1


a = -1
b = 1
N = 25
eq = SturmLiouville(p2, r2, f2, a, b, N, p_prime=p_prime2)
us_2, xs1_2 = eq.fem()

us2_2, xs_2 = eq.sem()


plt.plot(
    np.linspace(a, b, N + 1),
    [ana_solution2(x) for x in np.linspace(a, b, N + 1)],
    label="Analytical",
    alpha=0.5,
)
plt.plot(np.linspace(a, b, N + 1), us_2, label="FEM", linestyle=(0, (3, 10)))
plt.plot(xs_2, us2_2, label="SEM", linestyle=(1, (5, 10)))

plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
plt.grid()
plt.show()
plt.clf()
