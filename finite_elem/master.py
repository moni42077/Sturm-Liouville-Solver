import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
def phi1(i,xi):
    if i == 0:
        return 1 - xi
    elif i == 1:
        return xi
    else:
        return 0

def phi1_prime(i,xi):
    if i == 0:
        return -1
    elif i == 1:
        return 1
    else:
        return 0
    
    
def phi(i,x,xs):
    if i!=0 and xs[i-1] <= x <= xs[i]:
        return (x-xs[i-1])/(xs[i]-xs[i-1])
    elif i != len(xs) and xs[i] <= x <= xs[i+1]:
        return (xs[i+1]-x)/(xs[i+1]-xs[i])
    else:
        return 0 

def phi_prime(i,x,xs):
    if i!=0 and xs[i-1] <= x <= xs[i]:
        return 1/(xs[i]-xs[i-1])
    elif i != len(xs) and xs[i] <= x <= xs[i+1]:
        return 1/(xs[i+1]-xs[i])
    else:
        return 0 
    
    
N = 10
xs = np.linspace(-1,1,N+1)    
i = 1

A = np.zeros((N+1,N+1))
#Normal Phi
for i in range(N+1):
    if i == 0:
        f_ii = lambda x: phi_prime(i,x,xs)*phi_prime(i,x,xs)
        A[i][i] = quad(f_ii,xs[i],xs[i+1])[0]

        f_ii_p1 = lambda x: phi_prime(i,x,xs)*phi_prime(i+1,x,xs)
        A[i][i+1] = quad(f_ii_p1,xs[i],xs[i+1])[0]
    elif i == N:
        print("!!HELOO!!")
        f_ii = lambda x: phi_prime(i,x,xs)*phi_prime(i,x,xs)
        A[i][i] = quad(f_ii,xs[i-1],xs[i])[0]

        f_ii_m1 = lambda x: phi_prime(i,x,xs)*phi_prime(i-1,x,xs)
        A[i][i-1] = quad(f_ii_m1,xs[i-1],xs[i])[0]
    else:
        f_ii = lambda x: phi_prime(i,x,xs)*phi_prime(i,x,xs)
        A[i][i] = quad(f_ii,xs[i-1],xs[i+1])[0]

        f_ii_m1 = lambda x: phi_prime(i,x,xs)*phi_prime(i-1,x,xs)
        A[i][i+1] = quad(f_ii_m1,xs[i-1],xs[i])[0]

        f_ii_p1 = lambda x: phi_prime(i,x,xs)*phi_prime(i+1,x,xs)
        A[i][i-1] = quad(f_ii_p1,xs[i],xs[i+1])[0]

print(A)
print(f'1/h={1/(xs[7]-xs[6])}')

B = np.zeros((N+1,N+1))
for i in range(N+1):
    if i == N:
        h_i = xs[i]-xs[i-1]
    else:
        h_i = xs[i+1]-xs[i]
    J = h_i 
    x_transform = lambda xi: h_i * xi + xs[i]
    if i == 0:
        g_ii = lambda xi: phi1_prime(i,xi)*phi1_prime(i,xi)
        B[i][i] = quad(g_ii,0,1)[0] / J 
        
        g_ii_p1 = lambda xi: phi1_prime(i,xi)*phi1_prime(i+1,xi)
        B[i][i+1] = quad(g_ii,0,1)[0] / J 
    elif i == N:
        g_ii = lambda xi: phi1_prime(1,xi)*phi1_prime(1,xi)
        B[i][i] = quad(g_ii,0,1)[0] / J 
        
        g_ii_m1 = lambda xi: phi1_prime(1,xi)*phi1_prime(0,xi)
        B[i][i-1] = quad(g_ii,0,1)[0] / J 
    else:
        g_ii = lambda xi: phi1_prime(0,xi)*phi1_prime(0,xi)
        B[i][i] = quad(g_ii,0,1)[0] / J * 2
        
        g_ii_p1 = lambda xi: phi1_prime(0,xi)*phi1_prime(1,xi)
        B[i][i+1] = quad(g_ii,0,1)[0] / J 
        
        g_ii_m1 = lambda xi: phi1_prime(1,xi)*phi1_prime(0,xi)
        B[i][i-1] = quad(g_ii,0,1)[0] / J 
print(B)

i = 2
f = lambda x: phi(i,x,xs)*phi(i+1,x,xs)
print(quad(f,xs[i-1],xs[i+1])[0])