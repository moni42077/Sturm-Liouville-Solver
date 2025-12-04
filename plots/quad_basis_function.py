import matplotlib.pyplot as plt
import numpy as np


N = 3
xs = np.linspace(0,1,N)


def phi_stand_quad(i, xi):
    if i == 0:
        return 2*xi**2 - 3*xi + 1
    elif i == 1:
        return 4*xi*(1 - xi)
    elif i == 2:
        return 2*xi**2 - xi
    else:
        return 0
        
        
xs_dense = np.linspace(0,1,100)
plt.figure(figsize=(10,7))
plt.title(f'Master Quad Basis functions')
plt.xlabel(r'$\xi$')
plt.ylabel('y')
plt.grid()
plt.scatter(xs,[0]*N,label='Grid Points')
for i in range(N):
    plt.plot(xs_dense,[phi_stand_quad(i,xi) for xi in xs_dense],label=fr'$\varphi_{i}$')

#plt.legend(loc='upper right')
plt.legend()
#plt.show()


plt.savefig(f'quad_basis{N}.eps')