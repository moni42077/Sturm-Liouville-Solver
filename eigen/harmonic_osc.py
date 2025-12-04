import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from element import EigenvalueProblem


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


def f(x):
    return -1.0


N = 35
L = 5


a = 0
b = L

xs_cheb = np.array([np.cos(np.pi * j / (N)) for j in range(N + 1)])

if a == -1.0 and b == 1.0:
    xs = xs_cheb
else:
    xs = np.array([(b - a) * x / 2 + (a + b) / 2 for x in xs_cheb])

D = chebyshev_diff_matrix(xs, N)
A = D @ D
B = np.diag([f(x) for x in xs])

# lambda, eigenvalue
l, v = sp.linalg.eig(A[1:-1, 1:-1], B[1:-1, 1:-1])
l_ana = lambda x: (np.pi * x / L) ** 2
lam_ana_l = [l_ana(x) for x in range(10)]
for i in range(len(l) - 5, len(l)):
    if True:
        u = np.zeros(N + 1)
        u[1:-1] = v[:, i].real
        c = np.polyfit(xs, u, N)
        p = np.poly1d(c)
        x_dense = np.linspace(a, b, 100)
        y_dense = p(x_dense)

        # plt.plot(xs, u, label=f"Eigenvalue {l[i].real:.3f}")
        plt.plot(x_dense, y_dense, label=f"polyfit for {l[i].real:.3f}")
        sol = lambda x: np.sin(np.sqrt(l[i]) * x)
        # plt.plot(
        #     x_dense, [sol(x) for x in x_dense], label=f"Analytical for {l[i].real:.3f}"
        # )


plt.legend()
plt.xlabel("x")
plt.ylabel("u(x)")
plt.grid()
plt.show()

# FEM
p = lambda x: -1
r = lambda x: -1
p_prime = lambda x: 0

eq = EigenvalueProblem(p, r, a, b, N, p_prime)
l_fem, v_fem, x_fem = eq.fem_stand()
l_sem, v_sem, x_sem = eq.sem(ne=N,N=9)


l_ana = np.array([(np.pi * n / L) ** 2 for n in range(N)])
plt.clf()
plt.plot(np.sort(l_ana)[:30], label="Analytical")
plt.plot(np.sort(l.real)[:30], label="SM")
plt.plot(np.sort(l_fem.real)[:30], label="FEM")
plt.plot(np.sort(l_sem.real)[:30], label="SEM")

plt.title("Harmonic Oscillator eigenvalues")
plt.grid()
plt.xlabel("n")
plt.ylabel(r"$\lambda_n$")
plt.legend()
plt.show()

plt.clf()
# Sort and real
idl_sm = np.argsort(l)
l_sr = l[idl_sm].real
v_sr = v[idl_sm]
idl_fem = np.argsort(l_fem)
l_fem_sr = l_fem[idl_fem].real
v_fem_sr = v_fem[idl_fem]

idl_sem = np.argsort(l_sem)
l_sem_sr = l_sem[idl_sem].real
v_sem_sr = v_sem[idl_sem]

vals = [(np.pi * n / L) ** 2 for n in [1,2,3,4]]
eps = 0.1
for val in vals:
    # SM
    for i in range(len(l)):
        if np.abs(l[i].real - val) < eps:
            u = np.zeros(N + 1)
            u[1:-1] = v[:, i].real
            
            sol1 =  np.sin(np.sqrt(l[i].real) * xs)
            sol1 /= np.linalg.norm(sol1)
            if np.dot(u,sol1) < 0:
                u *= -1
                
                
            c = np.polyfit(xs, u, N)
            p = np.poly1d(c)
            xs = np.sort(xs)
            x_dense = np.linspace(a, b, 100)
            u_dense = p(x_dense)
            u_dense /= np.sqrt(np.trapezoid(u_dense**2,x_dense))
            print(f'Before normalization sol = {sol1} and norm = {np.trapezoid(sol1**2, xs)}') 
            sol1 /= np.sqrt(np.trapezoid(sol1**2,xs))
            print(f'After normalization sol = {sol1} and norm = {np.trapezoid(sol1**2, xs)}') 
            
            plt.plot(
                xs , sol1, label=f"Analytical for $\lambda$ = {l[i].real:.3f}",alpha=0.6
            )
            plt.plot(
                x_dense, u_dense, label=f"SM Polyfit for $\lambda$ = {l[i].real:.3f}",linestyle=(0, (5, 10))
            )
            u /=  np.sqrt(np.trapezoid(u**2,xs))
            sm_err = np.abs(u-sol1)
            e_sm = l[i].real
            break
    # FEM
    for i in range(len(l_fem)):
        if np.abs(l_fem[i].real - val) < eps:
            u_fem = np.zeros(N + 1)
            u_fem[1:-1] = v_fem[:, i].real
            u_fem /= np.sqrt(np.trapezoid(u_fem**2,x_fem))#np.linalg.norm(u_fem)
            sol1 = np.array([np.sin(np.sqrt(l_fem[i]) * x) for x in x_fem])
            sol1 /= np.linalg.norm(sol1)
            
            if np.dot(u_fem,sol1) < 0:
                u_fem *= -1
            plt.plot(x_fem, u_fem, label=f"FEM for $\lambda$ = {l_fem[i].real:.3f}",linestyle=(12, (5, 10)))
            sol1 /= np.sqrt(np.trapezoid(sol1**2,x_fem))
            fem_err = np.abs(u_fem-sol1)
            e_fem = l_fem[i].real
            break

#    SEM
    for i in range(len(l_sem)):
        if np.abs(l_sem[i].real - val) < eps:
            u_sem = np.zeros(len(l_sem) + 2)
            u_sem[1:-1] = v_sem[:, i].real
            u_sem /= np.sqrt(np.trapezoid(u_sem**2,x_sem))#np.linalg.norm(u_sem)
            sol1 = np.array([np.sin(np.sqrt(l_sem[i]) * x) for x in x_sem])
            sol1 /= np.linalg.norm(sol1) 
            if np.dot(u_sem,sol1) < 0:
                u_sem *= -1
            plt.plot(x_sem, u_sem, label=f"SEM for $\lambda$ = {l_sem[i].real:.3f}",linestyle=(3, (5, 10)))
            sol1 /= np.sqrt(np.trapezoid(sol1**2,x_sem))
            sem_err = np.abs(u_sem-sol1)
            e_sem = l_sem[i].real
            break
    
    #Plot
    plt.title("Eigenvector plot")
    plt.legend()
    plt.grid()
    plt.xlabel(r"$x$")
    plt.ylabel(r"$u(x)$")
    plt.show()
    plt.clf()
    
    plt.plot(xs,sm_err,label=f"SM for $\lambda$ = {e_sm:.3f}")
    plt.plot(x_fem,fem_err,label=f"FEM for $\lambda$ = {e_fem:.3f}")
    plt.plot(x_sem,sem_err,label=f"SEM for $\lambda$ = {e_sem:.3f}")
    plt.legend()
    plt.grid()
    plt.xlabel(r"$x$")
    plt.ylabel(r"$\bar{\epsilon}$")
    plt.yscale("log")
    plt.title("Eigenvector Error plot")
    plt.show()
    plt.clf()
    

