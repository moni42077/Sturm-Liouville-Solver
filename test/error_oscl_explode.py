import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from finite_elem.equation import SturmLiouville

# If True will save plots to a file
save = False



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
#N = 126
N = 1132
# a = -1
# b = 1
# N = 10
#polynomial degree
po = 1
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
us3, xs3 = eq.fem_stand()
us4, xs4 = eq.fem_stand_quad()


x_dense = np.linspace(a, b, 100)
# Error plot
# ? NOTE: You don't have to interpolate as we HAVE the function for the analytical solution
# us_sm = np.interp(xs,np.flip(xs2),np.flip(us2))
# plt.plot(xs,us_sm,label='SM interp',linestyle=(1, (5, 10)))

plt.title(f"ne = {N}")
plt.plot(x_dense, [ana_solution(x) for x in x_dense], label="Analytical", alpha=0.5)
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

# plt.figure(figsize=(10,6))
# sem_avg_err = []
# sem2_avg_err = []
# sem3_avg_err = []
# fem_mast_avg_err = []
# fem_mast_quad_avg_err = []
# sm_avg_err = []
# n_range = [i for i in range(3, 1000)]
# for n in n_range:
#     eq = SturmLiouville(p, r, f, a, b, n, p_prime=p_prime)
#     #us, xs = eq.sem(n,N=po)
#     #us_s2, xs_s2 = eq.sem(n,N=2*po)
#     #us_s3, xs_s3 = eq.sem(n,N=3*po+6)
#     #us2, xs2 = eq.sm()
#     us3, xs3 = eq.fem_stand()
#     us4, xs4 = eq.fem_stand_quad()
    
#     #sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
#     #sem2_err = [np.abs(ana_solution(x) - us_s2[i]) for i, x in enumerate(xs_s2)]
#     #sem3_err = [np.abs(ana_solution(x) - us_s3[i]) for i, x in enumerate(xs_s3)]
#    # sm_err = [np.abs(ana_solution(x) - us2[i]) for i, x in enumerate(xs2)]
#     fem_mast_err = [np.abs(ana_solution(x) - us3[i]) for i, x in enumerate(xs3)]
#     fem_mast_quad_err = [np.abs(ana_solution(x) - us4[i]) for i, x in enumerate(xs4)]
    
#     #sem_avg_err.append(np.mean(sem_err))
#     #sem2_avg_err.append(np.mean(sem2_err))
#     #sem3_avg_err.append(np.mean(sem3_err))
#     fem_mast_avg_err.append(np.mean(fem_mast_err))
#     fem_mast_quad_avg_err.append(np.mean(fem_mast_quad_err))
#     # if np.mean(fem_mast_quad_err) > 10:
#     #     print(f"For n = {n} error is {np.mean(fem_mast_quad_err)}")
#     #     plt.clf()
#     #     plt.plot(x_dense, [ana_solution(x) for x in x_dense], label="Analytical", alpha=0.5)
#     #     plt.plot(xs4,us4,label="Quad")
#     #     plt.legend()
#     #     plt.show()
#     #sm_avg_err.append(np.mean(sm_err))

# #plt.plot(n_range, sem_avg_err, label=f"SEM N={po}")
# #plt.plot(n_range, sem2_avg_err, label=f"SEM N={2*po}")
# #plt.plot(n_range, sem3_avg_err, label=f"SEM N={3*po+6}")
# plt.plot(n_range, fem_mast_avg_err, label="FEM Master")
# plt.plot(n_range, fem_mast_quad_avg_err, label="FEM Master Quad")
# #plt.plot(n_range, sm_avg_err, label="SM N = ne")
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
