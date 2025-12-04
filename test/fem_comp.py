import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from finite_elem.equation import SturmLiouville

# If True will save plots to a file
save = True


# def p(x):
#     return -11

# def p_prime(x):
#     return 0

# def r(x):
#     return 3*x

# def f(x):
#     #return np.sin(np.pi*x)
#     return -3*x**3*np.sin(np.pi*x) - 22*np.sin(np.pi*x)-44*x*np.pi*np.cos(np.pi*x)+11*x**2*np.pi**2*np.sin(np.pi*x)


# def ana_solution(x):
#     #return - np.sin(np.pi*x)/np.pi**2
#     # c1 = -1 * (np.exp(4) + np.exp(-4)) / 32
#     # c2 = -c1 - np.exp(4)/16
#     # return c1 + c2*x+np.exp(4*x)/16
#     return -x**2*np.sin(np.pi*x)

def p(x):
    return -3*x

def p_prime(x):
    return -3

def r(x):
    return 2

def f(x):
    return 2*x**3+25*x**2-28*x

def ana_solution(x):
    #u(-3)=u(2)=0
    return (x+3)*(x-2)**2


name = "eq2"
a = -3
b = 2
N = 10
#polynomial degree
fem_mast_avg_err = []
fem_avg_err = []
n_range = [i for i in range(3, 30)]
for n in n_range:
    eq = SturmLiouville(p, r, f, a, b, n, p_prime=p_prime)
    us, xs = eq.fem()
    us3,xs3 = eq.fem_stand()
    
    fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]
    fem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
    
    fem_mast_avg_err.append(np.mean(fem_mast_err))
    fem_avg_err.append(np.mean(fem_err))

plt.plot(n_range, fem_avg_err, label="FEM")
plt.plot(n_range, fem_mast_avg_err, label="FEM Master",linestyle=(1, (5, 10)))
plt.title("Log-Linear error comparison")
plt.xlabel(r"$ne$")
plt.ylabel(r"$\bar{\epsilon}$")
plt.legend()
plt.grid()
plt.yscale("log")
if save:
    plt.savefig(f"{name}_fem_comp.eps")
else:
    plt.show()
plt.clf()