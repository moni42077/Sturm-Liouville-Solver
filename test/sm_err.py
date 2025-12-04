
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from finite_elem.equation import SturmLiouville

# If True will save plots to a file
save = False


# def p(x):
#     return -1


# def p_prime(x):
#     return 0


# def r(x):
#     return 0


# def f(x):
#     return np.sin(np.pi * x)
#     # return 12*x**2*np.sin(x) + 8*x**3*np.cos(x) - x**4*np.sin(x)
#     # return -6*x**4+6*x**3+6*x**2-2


# def ana_solution(x):
#     return -np.sin(np.pi * x) / np.pi**2
#     # return -x**3+x**2+x-1
# #     return x**4*np.sin(x)


def p(x):
    return -11

def p_prime(x):
    return 0

def r(x):
    return 3*x

def f(x):
    #return np.sin(np.pi*x)
    return -3*x**3*np.sin(np.pi*x) - 22*np.sin(np.pi*x)-44*x*np.pi*np.cos(np.pi*x)+11*x**2*np.pi**2*np.sin(np.pi*x)

def ana_solution(x):
    #return - np.sin(np.pi*x)/np.pi**2
    # c1 = -1 * (np.exp(4) + np.exp(-4)) / 32
    # c2 = -c1 - np.exp(4)/16
    # return c1 + c2*x+np.exp(4*x)/16
    return -x**2*np.sin(np.pi*x)

# def p(x):
#     return -np.exp(-4*x)

# def p_prime(x):
#     return 4*np.exp(-4*x)

# def r(x):
#     return 4*np.exp(-4*x)

# def f(x):
#     c= -4 * np.e/(1+np.e**2)
#     return np.exp(-3*x)+c*np.exp(-4*x)

# def ana_solution(x):
#     c= -4 * np.e/(1+np.e**2)
#     return np.exp(x)-np.sinh(1)*np.exp(2*x)/np.sinh(2) + c/4





# def p(x):
#     return -1

# def p_prime(x):
#     return 0

# def r(x):
#     return 2*x

# def f(x):
#     return 2*np.cos(x)-x*np.sin(x)+2*x**2*np.sin(x)

# def ana_solution(x):
#     #-pi=pi =0
#     return x*np.sin(x)


name = "sem_x_sin"
a = -1
b = 1
N = 15
#polynomial degree
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
us, xs = eq.sm()
us2, xs2 = eq.sem(ne=1,N=N)


x_dense = np.linspace(a, b, 100)
# Error plot
# ? NOTE: You don't have to interpolate as we HAVE the function for the analytical solution
# us_sm = np.interp(xs,np.flip(xs2),np.flip(us2))
# plt.plot(xs,us_sm,label='SM interp',linestyle=(1, (5, 10)))

plt.title(f"N = {N}")
plt.plot(x_dense, [ana_solution(x) for x in x_dense], label="Analytical", alpha=0.5)
plt.plot(xs , us, label=f"SM N = {N}", linestyle=(0, (3, 10)))
plt.plot(xs2, us2, label=f"SEM ne = 1, N = {N}", linestyle=(1, (5, 10)))


plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.grid()
if save:
    plt.savefig(f"{name}.eps")
else:
    plt.show()
plt.clf()
#Error plot
sm_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
sm2_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]

plt.plot(xs, sm_err, label="SM error")
plt.plot(xs2, sm2_err, label="SM2 error")

plt.xlabel("x")
plt.ylabel(r"$\epsilon$")
plt.legend()
plt.grid()
if save:
    plt.savefig(f"{name}_error.eps")
else:
    plt.show()
plt.clf()

# Average error comparison
sm_avg_err = []
sm2_avg_err = []
sem_avg_err = []
n_range = [i for i in range(3, 20)]
for n in n_range:
    eq = SturmLiouville(p, r, f, a, b, n, p_prime=p_prime)
    us, xs = eq.sm()
    us3,xs3 = eq.sem(1,n)
    
    sm_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
    sem_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]
    
    sm_avg_err.append(np.mean(sm_err))
    sem_avg_err.append(np.mean(sem_err))

plt.plot(n_range, sm_avg_err, label=f"SM  error")
plt.plot(n_range, sem_avg_err, label=f"SEM ne = 1 error")
plt.title("Log-Linear error comparison")
plt.xlabel("N")
plt.ylabel(r"$\bar{\epsilon}$")
plt.legend()
plt.grid()
plt.yscale("log")
if save:
    plt.savefig(f"{name}_error_N.eps")
else:
    plt.show()
plt.clf()