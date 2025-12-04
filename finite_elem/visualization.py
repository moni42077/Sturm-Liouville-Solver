import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from equation import SturmLiouville


def p(x):
    return -1


def r(x):
    return 0


def f(x):
    return np.sin(np.pi * x)


def f_prime(x):
    return np.pi * np.cos(np.pi * x)


a = -1
b = 1
N = 5
eq = SturmLiouville(p, r, f, a, b, N)
xs = np.linspace(a, b, 100)

for i in range(5):
    plt.plot(xs, [eq.phi(i, x) for x in xs], label=f"phi_{i}")
plt.legend()
plt.show()

plt.clf()

plt.plot(xs, [f(x) for x in xs], label="f(x)")
for i in range(N):
    plt.plot(xs, [f(x) * eq.phi(i, x) for x in xs], label=rf"$f(x)\phi_{i}$")
plt.legend()
plt.title(r'Plot of $f(x)=\sin(\pi x)$ for $ne=5$')
plt.savefig('f(x)phi.eps')
#plt.show()

plt.clf()

plt.plot(xs, [f_prime(x) for x in xs], label="f_prime(x)")
for i in range(N):
    plt.plot(xs, [f_prime(x) * eq.phi(i, x) for x in xs], label=f"f(x)' * φ_{i}")
plt.legend()
plt.show()

plt.clf()


plt.plot(xs, [f(x) for x in xs], label="f(x)")
for i in range(N):
    plt.plot(xs, [f(x) * eq.phi_prime(i, x) for x in xs], label=f"f(x)φ'{i}")
plt.legend()
plt.show()

plt.clf()


# phi_i * phi_j

for i in range(2, 3):
    for j in range(N):
        plt.plot(
            xs, [eq.phi(i, x) * eq.phi(j, x) for x in xs], label=f"phi_{i} * φ_{j}"
        )

plt.legend()
plt.show()
plt.clf()

sum = 0
for i in range(N):
    integ = sp.integrate.quad(lambda x: f(x) * eq.phi(i, x), a, b)[0]
    sum += integ
    # print(integ)

print(f"Integral of f = {sp.integrate.quad(f,a,b)[0]}")
print(sum)
