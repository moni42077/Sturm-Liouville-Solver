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
def p(x):
    return -x

def p_prime(x):
    return -1

def r(x):
    return 1/x**3

def f(x):
    return np.cos(1/x)/x**2

def ana_solution(x):
    return np.sin(1/x)

name = "osc2"
a = 1/(16*np.pi)
b = 1/(4*np.pi)
N = 10
# a = -1
# b = 1
# N = 10
#polynomial degree
po = 1
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
us, xs = eq.sem(N,N=po)
us2, xs2 = eq.sm()
us3, xs3 = eq.fem_stand()
us4, xs4 = eq.fem_stand_quad()


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
plt.plot(xs4, us4, label="FEM Master Quad", linestyle=(3, (5, 10)))


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

# sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
# sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]
# fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]
# fem_mast_quad_err = [np.abs(ana_solution(x) - us4[i]) for i, x in enumerate(xs4)]

# plt.title(f"Error vs x for ne = {N}")
# plt.plot(xs, sem_err, label=f"SEM N={po}")
# plt.plot(xs2, sm_err, label="SM")
# plt.plot(xs3, fem_mast_err, label="FEM Master ")
# plt.plot(xs4, fem_mast_quad_err, label="FEM Master Quad")

# plt.xlabel("x")
# plt.ylabel(r"$\epsilon$")
# plt.legend()
# plt.grid()
# if save:
#     plt.savefig(f"{name}_error.eps")
# else:
#     plt.show()
# plt.clf()

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
    eq = SturmLiouville(p, r, f, a, b, n, p_prime=p_prime)
    us, xs = eq.sem(n,N=po)
    us_s2, xs_s2 = eq.sem(n,N=2*po)
    us_s3, xs_s3 = eq.sem(n,N=3*po+6)
    us2, xs2 = eq.sm()
    us3, xs3 = eq.fem_stand()
    us4, xs4 = eq.fem_stand_quad()
    
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

plt.plot(n_range, sem_avg_err, label=f"SEM N={po}")
plt.plot(n_range, sem2_avg_err, label=f"SEM N={2*po}")
plt.plot(n_range, sem3_avg_err, label=f"SEM N={3*po+6}")
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

plt.figure(figsize=(10,6))
sem_avg_err = []
sem2_avg_err = []
sem3_avg_err = []
fem_mast_avg_err = []
fem_mast_quad_avg_err = []
sm_avg_err = []
n_range = [i for i in range(3, 3000)]
for n in n_range:
    eq = SturmLiouville(p, r, f, a, b, n, p_prime=p_prime)
    #us, xs = eq.sem(n,N=po)
    #us_s2, xs_s2 = eq.sem(n,N=2*po)
    #us_s3, xs_s3 = eq.sem(n,N=3*po+6)
    #us2, xs2 = eq.sm()
    us3, xs3 = eq.fem_stand()
    us4, xs4 = eq.fem_stand_quad()
    
    #sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
    #sem2_err = [np.abs(ana_solution(x) - us_s2[i]) for i, x in enumerate(xs_s2)]
    #sem3_err = [np.abs(ana_solution(x) - us_s3[i]) for i, x in enumerate(xs_s3)]
   # sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]
    fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]
    fem_mast_quad_err = [np.abs(ana_solution(x) - us4[i]) for i, x in enumerate(xs4)]
    
    #sem_avg_err.append(np.mean(sem_err))
    #sem2_avg_err.append(np.mean(sem2_err))
    #sem3_avg_err.append(np.mean(sem3_err))
    fem_mast_avg_err.append(np.mean(fem_mast_err))
    fem_mast_quad_avg_err.append(np.mean(fem_mast_quad_err))
    if np.mean(fem_mast_quad_err) > 100:
        print(f"For n = {n} quad error is {np.mean(fem_mast_quad_err)}")
    
    if np.mean(fem_mast_err) > 100:
        print(f"For n = {n} lin error is {np.mean(fem_mast_err)}")
    #     plt.clf()
    #     plt.plot(x_dense, [ana_solution(x) for x in x_dense], label="Analytical", alpha=0.5)
    #     plt.plot(xs4,us4,label="Quad")
    #     plt.legend()
    #     plt.show()
    #sm_avg_err.append(np.mean(sm_err))

#plt.plot(n_range, sem_avg_err, label=f"SEM N={po}")
#plt.plot(n_range, sem2_avg_err, label=f"SEM N={2*po}")
#plt.plot(n_range, sem3_avg_err, label=f"SEM N={3*po+6}")
plt.plot(n_range, fem_mast_avg_err, label="FEM Master")
plt.plot(n_range, fem_mast_quad_avg_err, label="FEM Master Quad")
#plt.plot(n_range, sm_avg_err, label="SM N = ne")
plt.title("Log-Log error comparison")
plt.xlabel(r"$\log(ne)$")
plt.ylabel(r"$\bar{\epsilon}$")
plt.legend()
plt.grid()
plt.yscale("log")
plt.xscale("log")
if save:
    plt.savefig(f"{name}_error_logN.eps")
else:
    plt.show()
plt.clf()
