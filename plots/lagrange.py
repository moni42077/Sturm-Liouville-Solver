import matplotlib.pyplot as plt 
import numpy as np

def lagrange(i,x,xs):
    num = 1
    den = 1
    for k in range(len(xs)):
        if k != i:
            num *= (x-xs[k])
            den *= (xs[i]-xs[k])
    return num/den
            
N = 10       
xs = np.linspace(-1,1,N)

for i in range(N):
    plt.plot(xs,[lagrange(i,x,xs) for x in xs])

plt.show()