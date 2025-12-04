import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from finite_elem.equation import SturmLiouville

# If True will save plots to a file
save = True

# def p(x):
#     return -x**2/2

# def p_prime(x):
#     return -x

# def r(x):
#     return -3

# def f(x):
#     return -6*x

# def ana_solution(x):
#     #     a = -2
#     # A = -2
#     # b = 3
#     # B = 18
#     return x**2+3*x

def p(x):
    return -x**2/3

def p_prime(x):
    return -2*x/3

def r(x):
    return 1/3

def f(x):
    return 7*x**4-13*x**3+x

def ana_solution(x):
    #     a = -2
    # A = -2
    # b = 3
    # B = 18
    return x**4-3*x**3+x


def psi(x):
    return (B-A)*(x-a)/(b-a) + A





name = "non_hom_1"
a = -1
A = 3
b = 3
B = 3

# a = -2
# A = -2
# b = 3
# B = 18
N = 20
psi_prime =(B-A)/(b-a) 
g = lambda x:f(x) + p_prime(x)*psi_prime-r(x)*psi(x)
#polynomial degree
po = 1
eq = SturmLiouville(p, r, g, a, b, N, p_prime=p_prime)
vs, xs = eq.sem(N,N=po)
vs2, xs2 = eq.sm()
vs3, xs3 = eq.fem_stand()

us = [psi(x)+vs[i] for i,x in enumerate(xs)]
us2 = [psi(x)+vs2[i] for i,x in enumerate(xs2)]
us3 = [psi(x)+vs3[i] for i,x in enumerate(xs3)]
x_dense = np.linspace(a, b, 100)
# Error plot
# ? NOTE: You don't have to interpolate as we HAVE the function for the analytical solution
# us_sm = np.interp(xs,np.flip(xs2),np.flip(us2))
# plt.plot(xs,us_sm,label='SM interp',linestyle=(1, (5, 10)))

plt.title(f"ne = {N}")
plt.plot(x_dense, [ana_solution(x) for x in x_dense], label="Analytical", alpha=0.5)
plt.plot(xs , us, label=f"SEM N = {po}", linestyle=(0, (3, 10)))
plt.plot(xs2, us2, label="SM", linestyle=(1, (5, 10)))
plt.plot(xs3, us3, label="FEM Master Elem", linestyle=(2, (5, 10)))


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
sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]
fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]

plt.plot(xs, sem_err, label="SEM N={po}")
plt.plot(xs2, sm_err, label="SM")
plt.plot(xs3, fem_mast_err, label="FEM Master")

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
plt.figure(figsize=(10,6))
sem_avg_err = []
sem2_avg_err = []
sem3_avg_err = []
fem_mast_avg_err = []
fem_mast_quad_avg_err = []
sm_avg_err = []
n_range = [i for i in range(3, 30)]
for n in n_range:
    eq = SturmLiouville(p, r, g, a, b, n, p_prime=p_prime)
    
    vs, xs = eq.sem(n,N=po)
    vs_s2, xs_s2 = eq.sem(n,N=4)
    vs_s3, xs_s3 = eq.sem(n,N=5)
    vs2, xs2 = eq.sm()
    vs3, xs3 = eq.fem_stand()
    vs4, xs4 = eq.fem_stand_quad()
    
    us = [psi(x)+vs[i] for i,x in enumerate(xs)]
    us_s2 = [psi(x)+vs_s2[i] for i,x in enumerate(xs_s2)]
    us_s3 = [psi(x)+vs_s3[i] for i,x in enumerate(xs_s3)]
    us2 = [psi(x)+vs2[i] for i,x in enumerate(xs2)]
    us3 = [psi(x)+vs3[i] for i,x in enumerate(xs3)]
    us4 = [psi(x)+vs4[i] for i,x in enumerate(xs4)]
        
    
    
    sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
    sem2_err = [np.abs(ana_solution(x) - us_s2[i]) for i, x in enumerate(xs_s2)]
    sem3_err = [np.abs(ana_solution(x) - us_s3[i]) for i, x in enumerate(xs_s3)]
    sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]
    fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]
    fem_mast_quad_err = [np.abs(ana_solution(x) - us4[i]) for i, x in enumerate(xs4)]
    
    sem_avg_err.append(np.mean(sem_err))
    sem2_avg_err.append(np.mean(sem2_err))
    sem3_avg_err.append(np.mean(sem3_err))
    fem_mast_avg_err.append(np.mean(fem_mast_err))
    fem_mast_quad_avg_err.append(np.mean(fem_mast_quad_err))
    sm_avg_err.append(np.mean(sm_err))
    # vs, xs = eq.sem(n,N=po)
    # vs2, xs2 = eq.sm()
    # vs3, xs3 = eq.fem_stand()

    # us = [psi(x)+vs[i] for i,x in enumerate(xs)]
    # us2 = [psi(x)+vs2[i] for i,x in enumerate(xs2)]
    # us3 = [psi(x)+vs3[i] for i,x in enumerate(xs3)]
        
    # sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
    # sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]
    # fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]


plt.plot(n_range, sem_avg_err, label=f"SEM N={po}")
plt.plot(n_range, sem2_avg_err, label=f"SEM N={4}")
plt.plot(n_range, sem3_avg_err, label=f"SEM N={5}")
plt.plot(n_range, fem_mast_avg_err, label="FEM Master")
plt.plot(n_range, fem_mast_quad_avg_err, label="FEM Master Quad")
plt.plot(n_range, sm_avg_err, label="SM N=ne")
plt.title("Log-Linear error comparison")
plt.xlabel(r"$ne$")
plt.ylabel(r"$\bar{\epsilon}$")
plt.legend()
plt.grid()
plt.yscale("log")
if save:
    plt.savefig(f"{name}_error_N.eps")
else:
    plt.show()
plt.clf()