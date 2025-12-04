import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
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
    return 1.0


N = 42
a = -np.pi
b = np.pi

xs_cheb = np.array([np.cos(np.pi * j / (N)) for j in range(N + 1)])

if a == -1.0 and b == 1.0:
    xs = xs_cheb
else:
    xs = np.array([(b - a) * x / 2 + (a + b) / 2 for x in xs_cheb])

D = chebyshev_diff_matrix(xs, N)

# q = 0
# A = -1 * D @ D + 2 * q * np.diag([np.cos(2 * x) for x in xs])
# B = np.diag([f(x) for x in xs])

# # lambda, eigenvalue
# l, v = sp.linalg.eig(A[1:-1, 1:-1], B[1:-1, 1:-1])

# for i in range(10, 15):
#     if True:
#         u = np.zeros(N + 1)
#         u[1:-1] = v[:, i].real
#         c = np.polyfit(xs, u, N)
#         p = np.poly1d(c)

#         x_dense = np.linspace(a, b, 100)
#         y_dense = p(x_dense)

#         plt.plot(xs, u, label=f"Eigenvalue {l[i].real:.3f}")
#         # plt.plot(x_dense, y_dense, label=f"polyfit")


# plt.legend()
# plt.xlabel("x")
# plt.ylabel("u(x)")
# plt.grid()
# plt.show()


# # FEM
# p = lambda x: 1
# r = lambda x: 2 * q * np.cos(2 * x) - 1
# p_prime = lambda x: 0

# eq = EigenvalueProblem(p, r, a, b, N, p_prime)

# l_fem, v_fem, x_fem = eq.fem_stand()


# plt.clf()
# plt.plot(np.sort(l.real)[:-5], label="SM")
# plt.plot(np.sort(l_fem.real)[:-5], label="FEM")

# plt.title("Eigenvalues for Mathieu")
# plt.legend()
# plt.grid()
# plt.show()
sm_sol = []
qs = np.arange(0,15.2,0.2)
for q in qs:
    A = -1 * D @ D + 2 * q * np.diag([np.cos(2 * x) for x in xs])
    B = np.diag([f(x) for x in xs])

# lambda, eigenvalue
    l , v = sp.linalg.eig(A[1:-1, 1:-1], B[1:-1, 1:-1])
    sm_sol.append(np.sort(l.real)[:10])
    

sm_sol = np.array(sm_sol)

plt.figure(figsize=(6, 10)) 


for i in range(sm_sol.shape[1]):
    plt.plot(qs,sm_sol[:,i],label=rf'$\lambda_{{{i}}}$')



ax = plt.gca()
ax.yaxis.set_major_locator(MultipleLocator(4))
ax.xaxis.set_major_locator(MultipleLocator(5))


#Cool thing allows you to use all the values that the functions start at 
# all_eigs = sm_sol[0].flatten()
# unique_ticks = np.unique(np.round(all_eigs, 3))  
# plt.yticks(unique_ticks)



plt.title("Eigenvalues of Mathieu equation for different q")
plt.legend()
plt.xlabel("q")
plt.ylabel(r"$\lambda$")
plt.grid()
plt.show()


# # FEM
# p = lambda x: 1
# r = lambda x: 2 * q * np.cos(2 * x) - 1
# p_prime = lambda x: 0

# eq = EigenvalueProblem(p, r, a, b, N, p_prime)

# l_fem, v_fem, x_fem = eq.fem_stand()