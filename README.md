# Sturm-Liouville Equation Solver
## Introduction
This is the code for My MSc in Mathematics at the University of Leeds. It solved Sturm-Lioville equations of the form 

$$-\frac{d}{dx}(p(x)u'(x))+r(x)u(x)=f(x) \qquad u(a)=A\quad u(b)=B$$

by using the following three methods — FEM (Finite Element Method), SM (Spectral Method) and SEM (Spectral Element Method).
### Example Usage
 ```python 
def p(x):
   return -np.exp(-4*x)

def p_prime(x):
   return 4*np.exp(-4*x)

def r(x):
   return 4*np.exp(-4*x)

def f(x):
   c= -4 * np.e/(1+np.e**2)
   return np.exp(-3*x)+c*np.exp(-4*x)

a = 1
b = 1
# Number of elements
N = 10

#Polynomial Degree for SEM
po = 2

eq = SturmLiouville(p, r, f, a, b, N, p_prime=p_prime)
us, xs = eq.sem(N,N=po) #N elements and polynomials of degree po
us2, xs2 = eq.sm() #Assumes polynomial degree N
us3, xs3 = eq.fem_stand()
us4, xs4 = eq.fem_stand_quad()
 ``` 
## Main file 
A bit more cleaned up version of the main file for the BVP can be found in `main/equation.py`, the actual file used for calculations and plots is `finite_element/equation.py` - note it contains non-finished methods, only use the following functions given `eq=SturmLuioville()`
- `eq.fem()` - finite element method with linear basis
- `eq.fem_stand()` - finite element method using linear basis on master element
- `eq.fem_stand_quad()` - finite element method using quadratic basis on master element
- `eq.sm()` - spectral method using Chebyshev diffentiation matrix
- `eq.sm2()` - spectral method using Legendre polynoamials (it is equivalent to SEM with $ne=1$)
- `eq.sem()` - spectral element method that takes $ne$ - number of elements and $N$ - degree of polynomial

And there are some support functions that might be usefull
- `eq.chebyshev_diff_matrix()` - returns the Chebyshev differentiation matrix for given $xs$
- `eq.differentiate()` - calculates the derivative of a function using the Chebyshev differentiation matrix
- `eq.legendre_gll_nodes()` - returns the GLL interpolation nodes for given $N$
- `eq.gll_weights()` - returns the GLL weights for given $N$ and $xs$
- `eq.deriv_matrix_local()` - return the GLL derivative matrix for given $xs$

## Other files
The `test` folder contains a lot of the files that produce error plots such as `test/error.py` which produces all the erorr plots for BVP, note that in order to run those files one must do `python -m test.error` for the `error.py` file for example, this is so we can import the SturmLioville class from `equation.py`.

Due to the uniqueness of the Eigenvalue Problems every one has its own seperate file in `eigen`, which does all the calculations and typically produces the plots as well.

The rest of the folders are not really important but here is a quick rundown 
- `finite_elem` - the main file is there as well as other element method files
- `plots` - some helper files to make some more specific plots
