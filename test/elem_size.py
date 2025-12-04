import matplotlib.pyplot as plt
import numpy as np

from finite_elem.equation import SturmLiouville

# If True will save plots to a file
save = True

def p(x):
    return -1

def p_prime(x):
    return 0

def r(x):
    return 2*x

def f(x):
    return 2*np.cos(x)-x*np.sin(x) + 2*x**2*np.sin(x)

def ana_solution(x):
    #u(-3) = u(-2) = 0
    return x*np.sin(x)


# def p(x):
#     return -11

# def p_prime(x):
#     return 0

# def r(x):
#     return 3*x

# def f(x):
#     return -3*x**3*np.sin(np.pi*x) - 22*np.sin(np.pi*x)-44*x*np.pi*np.cos(np.pi*x)+11*x**2*np.pi**2*np.sin(np.pi*x)

# def ana_solution(x):
#     return -x**2*np.sin(np.pi*x)


name = "elem_size_xsinx"
a = -3*np.pi
b = 3*np.pi
N = 12
eq_lin = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime,xs_type='linear')
eq_cheb = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime,xs_type='cheb')


us_lin, xs_lin = eq_lin.fem()
us2_lin, xs2_lin = eq_lin.fem_stand()
us3_lin, xs3_lin = eq_lin.fem_stand_quad()

us_cheb, xs_cheb = eq_cheb.fem()
us2_cheb, xs2_cheb = eq_cheb.fem_stand()
us3_cheb, xs3_cheb = eq_cheb.fem_stand_quad()


fem_err_lin = [np.abs(ana_solution(x) - us_lin[i]) for i, x in enumerate(xs_lin)]
fem_mast_err_lin = [np.abs(ana_solution(x) - us2_lin[i]) for i, x in enumerate(xs2_lin)]
fem_mast_quad_err_lin = [np.abs(ana_solution(x) - us3_lin[i]) for i, x in enumerate(xs3_lin)]

fem_err_cheb = [np.abs(ana_solution(x) - us_cheb[i]) for i, x in enumerate(xs_cheb)]
fem_mast_err_cheb = [np.abs(ana_solution(x) - us2_cheb[i]) for i, x in enumerate(xs2_cheb)]
fem_mast_quad_err_cheb = [np.abs(ana_solution(x) - us3_cheb[i]) for i, x in enumerate(xs3_cheb)]


plt.plot(xs_lin,fem_err_lin,label='FEM homogenous points', linestyle=(0, (3, 10)))
plt.plot(xs_cheb,fem_err_cheb,label='FEM Cheb points', linestyle=(1, (5, 10)))
plt.xlabel("x")
plt.ylabel(r"$\epsilon$")
plt.legend()
plt.grid()

plt.show()
plt.clf()


plt.title('Error in Master Elements')
plt.plot(xs2_lin,fem_mast_err_lin,label='Homogenous points')
plt.plot(xs2_cheb,fem_mast_err_cheb,label='Chebyshev points')

plt.xlabel("x")
plt.ylabel(r"$\epsilon$")
plt.legend()
plt.grid()

if save:
    plt.savefig(f"{name}_cheb.eps")
else:
    plt.show()
plt.clf()


plt.title('Error in Master Quadratic Elements')
plt.plot(xs3_lin,fem_mast_quad_err_lin,label='Homogenous points')
plt.plot(xs3_cheb,fem_mast_quad_err_cheb,label='Chebyshev points')
plt.xlabel("x")
plt.ylabel(r"$\epsilon$")
plt.legend()
plt.grid()


if save:
    plt.savefig(f"{name}_cheb_quad.eps")
else:
    plt.show()
    
plt.clf()
