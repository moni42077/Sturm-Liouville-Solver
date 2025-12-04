import numpy as np
import scipy as sp

def f(x):
    return np.sin(x)


def int_lob(f,k):
    z = np.zeros(k+1)
    z[0] = -1.
    z[-1] = 1.

    w = np.zeros(k+1)
    w[0] = 2/(k*(k+1))
    w[-1] = w[0]

    if k > 1:
        zL,wL = sp.special.roots_legendre(k-1)
        for i in range(1,k):
            z[i] = zL[i-1]
            w[i] = wL[i-1]
            
    f_vals = [f(x) for x in z]
    result = 0
    for i in range(k+1):
        result += f_vals[i]*w[i]
    return result

f_int = lambda x: np.cos(x)

print(f'correct is {f_int(1) - f_int(-1)}')
print(int_lob(f,10))