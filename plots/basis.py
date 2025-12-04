import matplotlib.pyplot as plt
import numpy as np


N = 5
xs = np.linspace(0,1,N)

def phi_lin(i,x,xs):
    if i != 0 and xs[i-1] <= x <= xs[i]:
        return (x-xs[i-1])/(xs[i]-xs[i-1])
    elif i !=N and xs[i] <= x <= xs[i+1]: 
        return (xs[i+1]-x)/(xs[i+1]-xs[i])
    else:
        return 0

    


plt.figure(figsize=(10,7))
plt.title(f'Basis functions on a homogeneous interval with N={N}')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.scatter(xs,[0]*N,label='Grid Points')
for i in range(N):
    plt.plot(xs,[phi_lin(i,x,xs) for x in xs],label=fr'$\phi_{i}$')

plt.legend(loc='upper right')

#plt.show()

plt.savefig(f'basis{N}.eps')