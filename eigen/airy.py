import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from element import EigenvalueProblem
from matplotlib.ticker import MultipleLocator
from scipy.special import airy


def chebyshev_diff_matrix(xs, N):


    """Chebyshev differentiation matrix:
    Cite: Trefethen, Lloyd N. — Spectral Methods in MATLAB
    """
    # as shown in book p.54
    c = np.ones(N + 1)
    c[0] = 2
    c[-1] = 2
    c = c * (-1) ** np.arange(N + 1)  # c = [2,1,1,...,1,2]
    X = np.tile(xs, (N + 1, 1))  # [xs,xs,..,xs]
    dX = X - X.T  # diag 0
    D = (np.outer(c, 1 / c)) / (dX + np.eye(N + 1))  # eye diag 1s else 0
    D = D - np.diag(np.sum(D, axis=1))
    D = -1 * D
    return D


def b(x):
    return x


N = 28
xs_cheb = np.array([np.cos(np.pi * j / (N)) for j in range(N + 1)])

D = chebyshev_diff_matrix(xs_cheb, N)
#u''=l*x*u
A = D @ D
B = np.diag([b(x) for x in xs_cheb])

# lambda, eigenvalue
l, v = sp.linalg.eig(A[1:-1, 1:-1], B[1:-1, 1:-1])


#FEM
a = -1
b = 1
p = lambda x : -1
r = lambda x : -x
p_prime = lambda x : 0
eq = EigenvalueProblem(p,r,a,b,N,p_prime)
# l_fem,v_fem,x_fem = eq.fem_stand()

# for i in range(len(l_fem)):
#     if np.abs(l_fem[i].real - 500) <= 41:
#         u_fem = np.zeros(N + 1)
#         u_fem[1:-1] = v_fem[:, i].real
#         u_fem = np.flip(u_fem)
#         u_fem *= -1
#         u_fem /= np.sqrt(np.trapezoid(u_fem**2,x_fem))
#         sol2 = airy(l_fem[i].real**(1/3)*x_fem)[0]
#         sol2 /= np.sqrt(np.trapezoid(sol2**2,x_fem))
#         plt.plot(x_fem,sol2,label=f'Analytical Solution for $\lambda$={l_fem[i].real:.3f}',alpha=0.6)
#         plt.plot(x_fem,u_fem,label=f'FEM for $\lambda$={l_fem[i].real:.3f}',linestyle=(0,(5,10)))

#         plt.title(f"Airy Eigenvector plot for N = {N}")
#         plt.legend()
#         plt.grid()    
#         plt.show()
#         #plt.savefig(f"airy_fem_{N}.eps")
#         plt.clf()
        
#         #Err plot
#         err_lin_fem = np.abs(u_fem-sol2)
#         plt.plot(err_lin_fem)
#         plt.grid()
#         plt.title(f"Error for N = {N} and $\lambda$ = {l_fem[i].real:.3f}")
#         plt.xlabel(r"$n$")
#         plt.ylabel(r"$e_n$")
#         plt.show()
#         #plt.savefig(f"airy_fem_err_{N}.eps")
#         plt.clf()
       
       
       
l_sem,v_sem,x_sem = eq.sem(ne=N,N=9)

for i in range(len(l_sem)):
    if np.abs(l_sem[i].real - 500) <= 41:
        u_sem = np.zeros(9*N + 1)
        u_sem[1:-1] = v_sem[:, i].real
        u_sem = np.flip(u_sem)
        u_sem *= -1
        u_sem /= np.sqrt(np.trapezoid(u_sem**2,x_sem))
        sol2 = airy(l_sem[i].real**(1/3)*x_sem)[0]
        sol2 /= np.sqrt(np.trapezoid(sol2**2,x_sem))
        plt.plot(x_sem,sol2,label=f'Analytical Solution for $\lambda$={l_sem[i].real:.3f}',alpha=0.6)
        plt.plot(x_sem,u_sem,label=f'SEM N = 9 for $\lambda$={l_sem[i].real:.3f}',linestyle=(0,(5,10)))

        plt.title(f"Airy Eigenvector plot for ne = {N}")
        plt.legend()
        plt.grid()    
        plt.show()
        #plt.savefig(f"airy_fem_{N}.eps")
        plt.clf()
        
        #Err plot
        err_lin_sem = np.abs(u_sem-sol2)
        plt.plot(err_lin_sem)
        plt.grid()
        plt.title(f"SEM Error for ne = {N} and N = 9 and $\lambda$ = {l_sem[i].real:.3f}")
        plt.xlabel(r"$n$")
        plt.ylabel(r"$e_n$")
        plt.show()
        #plt.savefig(f"airy_fem_err_{N}.eps")
        plt.clf()       
       
        
"""Uncomment below for SM solution"""

# for i in range(len(l)):
#     if np.abs(l[i] - 500) <= 2:
#         u = np.zeros(N + 1)
#         u[1:-1] = v[:, i].real
#         u = np.flip(u)
#         u *= -1
#         xs = np.flip(xs_cheb)
#         u /= np.sqrt(np.trapezoid(u**2,xs))
#         c = np.polyfit(xs, u, N)
#         p = np.poly1d(c)

#         x_dense = np.linspace(-1, 1, 100)
#         y_dense = p(x_dense)

#         sol1 = airy(l[i].real**(1/3)*x_dense)[0]
#         sol1 /= np.sqrt(np.trapezoid(sol1**2,x_dense))#np.linalg.norm(sol1)
        
#         plt.plot(
#                 x_dense, sol1, label=f"Analytical for $\lambda$ = {l[i].real:.3f}",alpha=0.6
#             )
#         #plt.plot(xs, u, label=f"SM for $\lambda$ =  {l[i].real:.3f}",linestyle=(12,(5,10)))
#         plt.plot(x_dense, y_dense, label=f"SM polyfit for $\lambda$ = {l[i].real:.3f}",linestyle=(0,(5,10)))
        
#         ax = plt.gca()
#         ax.yaxis.set_major_locator(MultipleLocator(0.5))
#         ax.xaxis.set_major_locator(MultipleLocator(0.5))
#         plt.legend()
#         plt.title(f"Airy Eigenvector plot for N = {N}")
#         plt.xlabel(r"$x$")
#         # plt.xlim(-1,1)
#         # plt.ylim(-1,1)
#         plt.ylabel(r"$u(x)$")
#         plt.grid()
#         plt.show()
#         #plt.savefig(f"airy_eigen_{N}.eps")
#         plt.clf()
        
#         #Err at interpolation points
#         sol2 = airy(l[i].real**(1/3)*xs)[0]
#         sol2 /= np.sqrt(np.trapezoid(sol2**2,xs))#np.linalg.norm(sol1)
#         err_lin = np.abs(u-sol2)
#         #err_lin = np.sqrt(np.mean((u-sol2)**2))
#         plt.plot(err_lin)
        
#         plt.grid()
#         plt.title(f"Error for N = {N} and $\lambda$ = {l[i].real:.3f}")
#         plt.xlabel(r"$n$")
#         plt.ylabel(r"$e_n$")
#         #plt.savefig(f"airy_err_{N}.eps")
#         plt.show()
        

# plt.clf()

# plt.plot(np.sort(l.real)[5:30],label='SM')
# plt.plot(np.sort(l_fem.real)[5:30],label='FEM')

# # plt.scatter(np.sort(l[4:].real),np.sort(l_fem[4:].real),alpha=0.6)

# # real_vals = np.concatenate((l[4:].real,l_fem[4:].real))
# # l_min = real_vals.min()
# # l_max = real_vals.max()
# # plt.plot([l_min,l_max],[l_min,l_max],label='x = y')
# # plt.xlabel(r'SM $\lambda$')
# # plt.ylabel(r'FEM $\lambda$')

# plt.title("Airy Equation Eigenvalues")


# plt.legend()
# plt.grid()
# plt.show()