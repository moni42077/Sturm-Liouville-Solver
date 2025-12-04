import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from finite_elem.equation import SturmLiouville

# If True will save plots to a file
save = False


def p(x):
    return -1


def p_prime(x):
    return 0


def r(x):
    return 0


def f(x):
    return np.sin(np.pi * x)
    # return 12*x**2*np.sin(x) + 8*x**3*np.cos(x) - x**4*np.sin(x)
    # return -6*x**4+6*x**3+6*x**2-2


def ana_solution(x):
    return -np.sin(np.pi * x) / np.pi**2
    # return -x**3+x**2+x-1
#     return x**4*np.sin(x)

def f2(x):
    return x**4

def psi(j,p,x,x_pj):
    L_n = sp.special.legendre(N)
    L_n_prime = L_n.deriv()
    if x != x_pj:
        return (1-x**2)*L_n_prime(x)/(N*(N+1)*L_n(x_pj)*(x-x_pj))
    else:
        return 1
    

a = -1
b = 1
N = 10
eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
xs = eq.legendre_gll_nodes()
ws = eq.gll_weights(xs)
ne = 5
int1 = 0 
int2 = 0

element_nodes = np.linspace(-1,1,ne+1)
l = np.diff(element_nodes)
for j in range(ne):
    J = l[j]/2
    intj1 = 0 
    intj2 = 0 
    
    xj1 = np.array([ws[j]+l[j]*(xs[i]+1)/2 for i in range(N+1)])
    xj2 = np.array([element_nodes[j]+l[j]*(xs[i]+1)/2 for i in range(N+1)])

    wj = np.array([l[j]*ws[i]/2 for i in range(N+1)])
    
    f_vals1 = np.array([f2(x) for x in xj1])
    f_vals2 = np.array([f2(x) for x in xj2])
    
    int1 += J*np.sum(wj*f_vals1)
    int2 += np.sum(wj*f_vals2)
    
    print(f'Intj1 is {int1} and inj2 is {int2} where quad is {sp.integrate.quad(f2,element_nodes[j],element_nodes[j+1])[0]}')
    
quad = 0 
for i in range(N):
    quad += f2(xs[i])*ws[i]

print(f'Int 1 is {int1} and int2 is {int2} and quadrature is {quad} and quad is {sp.integrate.quad(f2,-1,1)[0]}')

#f_product = lambda x: f2(x)*sp.special.eval_legendre(N,x)
#psi_vals2 = np.array([psi(j,p,x,xj1[p]) for p,x in enumerate(xj2)])