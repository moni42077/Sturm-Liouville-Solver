import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import legendre, roots_legendre


# Cheb polynoimial
def cheb_polyn(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return 2 * x * cheb_polyn(n - 1, x) - cheb_polyn(n - 2, x)


def cardinal(j, x, xs):
    p = np.ones(N + 1)
    p[0] = 2
    p[N] = 2
    return (
        2
        / (N * p[j])
        * np.sum([cheb_polyn(m, xs[j]) * cheb_polyn(m, x) / p[m] for m in range(N + 1)])
    )

# def cardinal_cheb_prime(i,j,x,xs):
    #d/dxi(Cj)
# cheb lobbato grid
N = 10
xs = [np.cos(np.pi * i / N) for i in range(N + 1)]
# xs.reverse()
j = 0
xs_dense = np.linspace(-1, 1, 100)

fig, axs = plt.subplots(3)

fig.suptitle("Chebyshev polynomial cardinal functions plot N = 10")
axs[0].plot(xs, [0] * len(xs), "o")
axs[1].plot(xs, [0] * len(xs), "o")
axs[2].plot(xs, [0] * len(xs), "o")
axs[0].plot(xs_dense, [cardinal(3, x, xs) for x in xs_dense], label="j=3")
axs[1].plot(xs_dense, [cardinal(5, x, xs) for x in xs_dense], label="j=5")
axs[2].plot(xs_dense, [cardinal(7, x, xs) for x in xs_dense], label="j=7")
axs[0].legend()
axs[1].legend()
axs[2].legend()
plt.show()


def gegenbauer_polyn(m, n, x):
    if n == 0:
        return 1
    elif n == 1:
        return 2 * m * x
    else:
        return (
            2 * (n - 1 + m) * x * gegenbauer_polyn(m, n - 1, x)
            - (n + 2 * m - 2) * gegenbauer_polyn(m, n - 2, x)
        ) / n


def legendre_polyn(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return (
            (2 * n - 1) * x * legendre_polyn(n - 1, x)
            - (n - 1) * legendre_polyn(n - 2, x)
        ) / n


def legendre_first_der(n, x):
    return gegenbauer_polyn(1.5, n - 1, x)


def cardinal_legendre(j, x, xs):
    return (
        -(1 - x**2)
        * legendre_first_der(N, x)
        / (N * (N + 1) * legendre_polyn(N, xs[j]) * (x - xs[j]))
    )


legendre_der_N = lambda x: legendre_first_der(N, x)

print(brentq(legendre_der_N, -1, 1))

fig, axs = plt.subplots(3)
P = legendre(N)
P_der = P.deriv()
f = P_der
y_dense = f(xs_dense)
sign_changes = np.where(np.diff(np.sign(y_dense)))[0]
roots = []
for i in sign_changes:
    x0, x1 = xs_dense[i], xs_dense[i + 1]
    try:
        root = brentq(f, x0, x1)
        roots.append(root)
    except:
        continue

xs = roots
print(f'Roots are {roots}')
print(f'Sympy roots are {roots_legendre(N-1)[0]}')
xs.insert(0, -1)
xs.append(1)
fig.suptitle("Legendre polynomial cardinal functions plot N = 10")
axs[0].plot(xs, [0] * len(xs), "o")
axs[1].plot(xs, [0] * len(xs), "o")
axs[2].plot(xs, [0] * len(xs), "o")
axs[0].plot(xs_dense, [cardinal_legendre(3, x, xs) for x in xs_dense], label="j=3")
axs[1].plot(xs_dense, [cardinal_legendre(5, x, xs) for x in xs_dense], label="j=5")
axs[2].plot(xs_dense, [cardinal_legendre(7, x, xs) for x in xs_dense], label="j=7")
axs[0].legend()
axs[1].legend()
axs[2].legend()
plt.show()
