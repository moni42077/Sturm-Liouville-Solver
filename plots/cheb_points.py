import matplotlib.pyplot as plt
import numpy as np


N = 20
xs_cheb = np.array([np.cos(np.pi * j / (N)) for j in range(N)])


plt.title(f'Chebyshev points for N={N}')
plt.xlabel('x')
plt.ylabel('y')
plt.scatter(xs_cheb,[0]*N)
plt.grid()
plt.show()

#plt.savefig(f'cheb{N}.eps')