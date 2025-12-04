import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from finite_elem.equation import SturmLiouville

# If True will save plots to a file
save = True


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
#     return 0

# def f(x):
#     return 2*np.cos(x)-x*np.sin(x)

# def ana_solution(x):
#     #-pi=pi =0
#     return x*np.sin(x)



# def p(x):
#     return -3*x

# def p_prime(x):
#     return -3

# def r(x):
#     return 2

# def f(x):
#     return 2*x**3+25*x**2-28*x

# def ana_solution(x):
#     #u(-3)=u(2)=0
#     return (x+3)*(x-2)**2




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


#Oscilating function
# def p(x):
#     return -x

# def p_prime(x):
#     return -1

# def r(x):
#     return 1/x**3

# def f(x):
#     return np.cos(1/x)/x**2

# def ana_solution(x):
#     return np.sin(1/x)

name = "r1_c"
# a = 1/(16*np.pi)
# b = 1/(4*np.pi)
# N = 10
a = -1
b = 1
#polynomial degree
po = 1



# Average error comparison
plt.figure(figsize=(10,6))
sem_avg_err = []
sem2_avg_err = []
sem3_avg_err = []
fem_mast_avg_err = []
fem_mast_quad_avg_err = []
sm_avg_err = []
n_range = [i for i in range(3, 200)]
for n in n_range:
    eq = SturmLiouville(p, r, f, a, b, n, p_prime=p_prime)
    us, xs = eq.sem(n,N=po)
    us3, xs3 = eq.fem_stand()
    
    sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
    
    fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]
    
    sem_avg_err.append(np.mean(sem_err))
    
    fem_mast_avg_err.append(np.mean(fem_mast_err))

    if n < 101:
        us_s2, xs_s2 = eq.sem(n,N=2*po)
        us4, xs4 = eq.fem_stand_quad()
        
        fem_mast_quad_err = [np.abs(ana_solution(x) - us4[i]) for i, x in enumerate(xs4)]
        fem_mast_quad_avg_err.append(np.mean(fem_mast_quad_err)) 
        sem2_err = [np.abs(ana_solution(x) - us_s2[i]) for i, x in enumerate(xs_s2)]
        sem2_avg_err.append(np.mean(sem2_err))
    
    if n < 26: 
        us_s3, xs_s3 = eq.sem(n,N=3*po+6)
        us2, xs2 = eq.sm()
        
        sem3_err = [np.abs(ana_solution(x) - us_s3[i]) for i, x in enumerate(xs_s3)]
        sem3_avg_err.append(np.mean(sem3_err))
        sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]
        sm_avg_err.append(np.mean(sm_err))

ne1 = po
ne2 = 2*po
ne3 = 3*po + 6
c_sem = [3*ne1*n-5 for n in n_range]
c_sem2 = [3*ne2*n-5 for n in n_range]
c_sem3 = [3*ne3*n-5 for n in n_range]

c_fem = [3*n-5 for n in n_range]
c_fem2 = [6*n-5 for n in n_range]

c_sm = [(n-1)**2 for n in n_range]

plt.plot(c_sem, sem_avg_err, label=f"SEM N={po}")
plt.plot(c_sem2[:len(sem2_avg_err)], sem2_avg_err, label=f"SEM N={2*po}")
plt.plot(c_sem3[:len(sem3_avg_err)], sem3_avg_err, label=f"SEM N={3*po+6}")
plt.plot(c_fem, fem_mast_avg_err, label="FEM Master")
plt.plot(c_fem2[:len(fem_mast_quad_avg_err)], fem_mast_quad_avg_err, label="FEM Master Quad")
plt.plot(c_sm[:len(sm_avg_err)], sm_avg_err, label="SM N=ne")
plt.title("Log-Linear error comparison")
plt.xlabel(r"$c$")
plt.ylabel(r"$\bar{\epsilon}$")
plt.legend()
plt.grid()
plt.yscale("log")
if save:
    plt.savefig(f"{name}_error_N.eps")
else:
    plt.show()
# plt.clf()

# plt.figure(figsize=(10,6))
# sem_avg_err = []
# sem2_avg_err = []
# sem3_avg_err = []
# fem_mast_avg_err = []
# fem_mast_quad_avg_err = []
# sm_avg_err = []
# n_range = [i for i in range(3, 100)]
# for n in n_range:
#     eq = SturmLiouville(p, r, f, a, b, n, p_prime=p_prime)
#     us, xs = eq.sem(n,N=po)
#     us_s2, xs_s2 = eq.sem(n,N=2*po)
#     us_s3, xs_s3 = eq.sem(n,N=3*po+6)
#     us2, xs2 = eq.sm()
#     us3, xs3 = eq.fem_stand()
#     us4, xs4 = eq.fem_stand_quad()
    
#     sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
#     sem2_err = [np.abs(ana_solution(x) - us_s2[i]) for i, x in enumerate(xs_s2)]
#     sem3_err = [np.abs(ana_solution(x) - us_s3[i]) for i, x in enumerate(xs_s3)]
#     sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]
#     fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]
#     fem_mast_quad_err = [np.abs(ana_solution(x) - us4[i]) for i, x in enumerate(xs4)]
    
#     sem_avg_err.append(np.mean(sem_err))
#     sem2_avg_err.append(np.mean(sem2_err))
#     sem3_avg_err.append(np.mean(sem3_err))
#     fem_mast_avg_err.append(np.mean(fem_mast_err))
#     fem_mast_quad_avg_err.append(np.mean(fem_mast_quad_err))
#     sm_avg_err.append(np.mean(sm_err))

# plt.plot(n_range, sem_avg_err, label=f"SEM N={po}")
# plt.plot(n_range, sem2_avg_err, label=f"SEM N={2*po}")
# plt.plot(n_range, sem3_avg_err, label=f"SEM N={3*po+6}")
# plt.plot(n_range, fem_mast_avg_err, label="FEM Master")
# plt.plot(n_range, fem_mast_quad_avg_err, label="FEM Master Quad")
# plt.plot(n_range, sm_avg_err, label="SM N = ne")
# plt.title("Log-Log error comparison")
# plt.xlabel(r"$\log(ne)$")
# plt.ylabel(r"$\bar{\epsilon}$")
# plt.legend()
# plt.grid()
# plt.yscale("log")
# plt.xscale("log")
# if save:
#     plt.savefig(f"{name}_error_logN.eps")
# else:
#     plt.show()
# plt.clf()
