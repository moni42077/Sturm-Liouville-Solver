import matplotlib.pyplot as plt
import numpy as np


N = 2
xs = np.linspace(0,1,N)

def phi_m(i,xi):
    if i == 0:
        return 1 - xi
    elif i == 1:
        return xi
    else:
        return 0

    


plt.figure(figsize=(10,7))
plt.title(f'Master Linear Basis functions')
plt.xlabel(r'$\xi$')
plt.ylabel('y')
plt.grid()
plt.scatter(xs,[0]*N,label='Grid Points')
for i in range(N):
    plt.plot(xs,[phi_m(i,xi) for xi in [0,1]],label=fr'$\phi_{i}$')

plt.legend(loc='upper right')

plt.show()

#plt.savefig(f'basis{N}.eps')