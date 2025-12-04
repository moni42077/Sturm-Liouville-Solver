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


a = -1.0
b = 1.0

xs_cheb = np.array([np.cos(np.pi * j / (N)) for j in range(N + 1)])

if a == -1.0 and b == 1.0:
    xs = xs_cheb
else:
    xs = np.array([(b - a) * x / 2 + (a + b) / 2 for x in xs_cheb])

D = chebyshev_diff_matrix(xs, N)

g = np.diag([1 - x**2 for x in xs])
A = D @ g @ D
B = np.diag([f(x) for x in xs])

# lambda, eigenvalue
l, v = sp.linalg.eig(A[1:-1, 1:-1], B[1:-1, 1:-1])
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
        # plt.plot(
        #     x_dense, [sol(x) for x in x_dense], label=f"Analytical for {l[i].real:.3f}"
        # )

plt.legend()
plt.xlabel("x")
plt.ylabel("u(x)")
plt.grid()
plt.show()
plt.clf()

# FEM
p = lambda x: 1 - x**2
r = lambda x: 1
p_prime = lambda x: -2 * x

eq = EigenvalueProblem(p, r, a, b, N, p_prime)

l_fem, v_fem, x_fem = eq.fem_stand()
l_sem, v_sem, x_sem = eq.sem(ne=N, N=9)

l_ana = np.array([n * (n + 1) for n in range(N)])

plt.clf()
plt.plot(np.sort(l.real)[:30], label="SM")
plt.plot(np.sort(l_fem.real)[:30], label="FEM")
plt.plot(np.sort(l_sem.real)[:30], label="SEM")
plt.plot(np.sort(l_ana)[:30], label="Analytical")

plt.title("Eigenvalues for Legendre")
plt.legend()
plt.grid()
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

vals = [2, 6]
eps = 0.5
eps = 0.1
for val in vals:
    # SM
    for i in range(len(l)):
        if np.abs(l[i].real - val) < eps:
            u = np.zeros(N + 1)
            u[1:-1] = v[:, i].real
            u /= np.linalg.norm(u)
            
            sol1 = np.array([np.sin(np.sqrt(l[i]) * x) for x in xs])
            sol1 /= np.linalg.norm(sol1)
            if np.dot(u,sol1) < 0:
                u *= -1
                
                
            c = np.polyfit(xs, u, N)
            p = np.poly1d(c)
            print(u)
            x_dense = np.linspace(a, b, 100)
            u_dense = p(x_dense)

            plt.plot(
                x_dense, u_dense, label=f"SM Polyfit for $\lambda$ = {l[i].real:.3f}"
            )
            
            plt.plot(
                xs , sol1, label=f"Analytical for {l[i].real:.3f}"
            )
            break
    # FEM
    for i in range(len(l_fem)):
        if np.abs(l_fem[i].real - val) < eps:
            u_fem = np.zeros(N + 1)
            u_fem[1:-1] = v_fem[:, i].real
            u_fem /= np.linalg.norm(u_fem)
            sol1 = np.array([np.sin(np.sqrt(l_fem[i]) * x) for x in x_fem])
            sol1 /= np.linalg.norm(sol1)
            
            if np.dot(u_fem,sol1) < 0:
                u_fem *= -1
            plt.plot(x_fem, u_fem, label=f"FEM for $\lambda$ = {l_fem[i].real:.3f}")
            break
    # for i in range(len(l_fem_sr)):
    #     if np.abs(l_fem_sr[i] - val) < eps:
    #         u_fem = np.zeros(N + 1)
    #         u_fem[1:-1] = v_fem_sr[:, i].real
    #         plt.plot(x_fem, u_fem, label=f"FEM for $\lambda$ = {l_fem_sr[i]:.3f}")
    #         break

#    SEM
    for i in range(len(l_sem)):
        if np.abs(l_sem[i].real - val) < eps:
            u_sem = np.zeros(len(l_sem) + 2)
            u_sem[1:-1] = v_sem[:, i].real
            u_sem /= np.linalg.norm(u_sem)
            sol1 = np.array([np.sin(np.sqrt(l_sem[i]) * x) for x in x_sem])
            sol1 /= np.linalg.norm(sol1) 
            if np.dot(u_sem,sol1) < 0:
                u_sem *= -1
            plt.plot(x_sem, u_sem, label=f"SEM for $\lambda$ = {l_sem[i].real:.3f}")
    
    #Analytical
    plt.legend()
    plt.show()
    plt.clf()
    # for i in range(len(l_sem_sr)):
    #     if np.abs(l_sem_sr[i] - val) < eps:
    #         u_sem = np.zeros(len(l_sem_sr) + 2)
    #         u_sem[1:-1] = v_sem_sr[:, i].real
    #         plt.plot(x_sem, u_sem, label=f"SEM for $\lambda$ = {l_sem_sr[i]:.3f}")

    # Analytical
    # val = 2 -> n = 1


plt.legend()
plt.xlabel("x")
plt.ylabel("u(x)")
plt.grid()
plt.show()


# print(f"SM is \n {l_sr} \n and SEM is {l_sem_sr}")
