import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from finite_elem.equation import SturmLiouville

# If True will save plots to a file
save = True

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


a = -1
b = 1
sem_errs = []
nes = []
Ns = []
for ne in range(10,20):
    for N in range(3,10):
        eq = SturmLiouville(p, r, f, a, b, ne, p_prime=p_prime)
        us, xs = eq.sem4(ne,N=N)
        sem_err = [np.abs(ana_solution(x) - us[i]) for i, x in enumerate(xs)]
        sem_errs.append(np.mean(sem_err))
        nes.append(ne)
        Ns.append(N)
        

#3D plot
# ax = plt.figure().add_subplot(projection='3d')

# ax.scatter(nes,Ns,sem_errs)
# ax.set_xlabel('ne')
# ax.set_ylabel('N')
# ax.set_zlabel(r'$\epsilon$')

#2D colormap
scatter = plt.scatter(nes,Ns,c=sem_errs,cmap='viridis')
plt.colorbar(scatter,label=r'$\epsilon$')

plt.xlabel('ne')
plt.ylabel('N')

plt.show()