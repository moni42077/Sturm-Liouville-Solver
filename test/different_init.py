import matplotlib.pyplot as plt
import numpy as np

from finite_elem.equation import SturmLiouville


def p(x):
    return 1


def p_prime(x):
    return 0


def r(x):
    return 0


def f(x):
    return np.cos(np.pi * x) * np.pi**2


def ana_solution(x):
    return np.cos(np.pi * x)


# u(a) = u(b) = 0
a = -1.5
b = 1.5
N = 25
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
us, xs = eq.fem()
us2, xs2 = eq.sem()


plt.plot(xs, [ana_solution(x) for x in xs], label="Analytical", alpha=0.5)
plt.plot(xs, us, label="FEM", linestyle=(0, (3, 10)))
plt.plot(xs2, us2, label="SM", linestyle=(1, (5, 10)))


plt.title(r"$-u(x)'' =\pi^2 \cos(\pi x))$ where $u(-1.5) = u(1.5) = 0$")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
plt.grid()
plt.savefig("cos.eps")  # plt.show()
plt.clf()

# Error plot
# Notice xs != xs2 so we want to transform xs2,us2 -> xs,us_sm
fem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]

plt.plot(xs, fem_err, label="FEM error")
plt.plot(xs, sm_err, label="SM error")

plt.title("|u_analytical(x)-u(x)|")
plt.xlabel("x")
plt.ylabel("Error")
plt.legend()
plt.grid()
plt.savefig("cos_error.eps")  # plt.show()
plt.clf()


def p(x):
    return -1


def p_prime(x):
    return 0


def r(x):
    return 0


def f(x):
    return 2 * (1 - 4 * x**2 * np.tanh(x**2)) / (np.cosh(x**2) ** 2)


def ana_solution(x):
    return np.tanh(x**2) - np.tanh(4)


# u(a) = u(b) = 0
a = -2
b = 2
N = 25
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
us, xs = eq.fem()
us2, xs2 = eq.sem()


plt.plot(xs, [ana_solution(x) for x in xs], label="Analytical", alpha=0.5)
plt.plot(xs, us, label="FEM", linestyle=(0, (3, 10)))
plt.plot(xs2, us2, label="SM", linestyle=(1, (5, 10)))


plt.title(r"$u(x)'' = 2(1-4x^2\tanh(x^2))\text{sech}^2(x^2)$ where $u(-2) = u(2) = 0$")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
plt.grid()
plt.savefig("tanh.eps")  # plt.show()
plt.clf()

# Error plot
# Notice xs != xs2 so we want to transform xs2,us2 -> xs,us_sm
fem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]

plt.plot(xs, fem_err, label="FEM error")
plt.plot(xs, sm_err, label="SM error")

plt.title("|u_analytical(x)-u(x)|")
plt.xlabel("x")
plt.ylabel("Error")
plt.legend()
plt.grid()
plt.savefig("tanh_error.eps")  # plt.show()
plt.clf()

# Non symetric


def p(x):
    return -1


def p_prime(x):
    return 0


def r(x):
    return 0


def f(x):
    return 6 * x - 5


def ana_solution(x):
    return x**3 - 2.5 * x**2 + 0.5 * x + 1


# u(a) = u(b) = 0
a = -0.5
b = 2
N = 25
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
us, xs = eq.fem()
us2, xs2 = eq.sem()


plt.plot(xs, [ana_solution(x) for x in xs], label="Analytical", alpha=0.5)
plt.plot(xs, us, label="FEM", linestyle=(0, (3, 10)))
plt.plot(xs2, us2, label="SM", linestyle=(1, (5, 10)))


plt.title(r"$u(x)'' = 6x-5$ where $u(-0.5) = u(2) = 0$")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
plt.grid()
plt.savefig("polyn2.eps")  # plt.show()
plt.clf()

# Error plot
# Notice xs != xs2 so we want to transform xs2,us2 -> xs,us_sm
fem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]

plt.plot(xs, fem_err, label="FEM error")
plt.plot(xs, sm_err, label="SM error")

plt.title("|u_analytical(x)-u(x)|")
plt.xlabel("x")
plt.ylabel("Error")
plt.legend()
plt.grid()
plt.savefig("polyn2_error.eps")  # plt.show()
plt.clf()
