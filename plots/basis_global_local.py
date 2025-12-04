import matplotlib.pyplot as plt
import numpy as np

def phi(i,xi):
    if i == 0:
        return 1-xi
    elif i == 1:
        return xi
    else:
        return 0

def x_inv_transform(i,x,xs):
    h_i = xs[i+1] - xs[i]
    return (x-xs[i])/h_i

def phi_transform(i,x,xs):
    if i == 0:
        if xs[0] <= x <= xs[1]:
            return phi(0,x_inv_transform(0,x,xs))
        else:
            return 0 
    elif i == N:
        if xs[N-1] <= x <= xs[N]:
            return phi(1,x_inv_transform(N-1,x,xs))
    else:
        if xs[i-1] <= x <= xs[i]:
            return phi(1,x_inv_transform(i-1,x,xs))
        elif xs[i] <= x <= xs[i+1]:
            return phi(0,x_inv_transform(i,x,xs))
        else:
            return 0
        
N = 3
a = 0
b = 1
#global nodes
xs = np.linspace(a,b,N+1)


plt.scatter(xs,[0]*(N+1),color='black',label='Points')

xs_dense = np.linspace(a,b,100)

for i in range(N+1):
    plt.plot(xs_dense,[phi_transform(i,x,xs) for x in xs_dense],label=rf'$\phi_{i}$')



plt.xticks(ticks=xs, labels=[r"$x_0$", r"$x_1$", r"$x_2$", r"$x_3$"])
plt.legend(loc='upper right')
#plt.savefig('basis_global_local.eps')
plt.show()

