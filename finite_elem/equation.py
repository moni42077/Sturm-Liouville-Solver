from types import FunctionType
import numpy as np
import scipy as sp
from scipy.fft import fft, ifft
from scipy.linalg import solve
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve


class SturmLiouville:
    """Sturm-Liouville equation of the form -d/dx(p(x)du/dx) + r(x)u = f(x)"""

    def __init__(
        self,
        p: FunctionType,
        r: FunctionType,
        f: FunctionType,
        a: int,
        b: int,
        N: int,
        p_prime=None,
        xs_type = 'Linear'
    ):
        """Initialize Equation  -d/dx(p(x)du/dx) + r(x)u = f(x)

        Args:
            p (FunctionType): p(x)
            r (FunctionType): r(x)
            f (FunctionType): f(x)
            a (int): u(a) = 0
            b (int): u(b) = 0
            N (int): number of points,
            p_prime (FunctionType) : p'(x) if available, to increase accuracy,
            xs_type (string) : spacing between the nodes for FEM, by default they are equispaced
        """
        self.p = p
        self.r = r
        self.f = f
        # make sure a is the smaller one
        if a <= b:
            self.a = a
            self.b = b
        else:
            self.a = b
            self.b = a

        self.N = N
        if xs_type == 'Linear' or xs_type =='linear':
            self.x = np.linspace(a, b, N + 1)
            
            #uncommet below if you want more concentrated points at the start (used in oscillating function example)
            # half = N // 2
            # quarter_point = (3*self.a+self.b)/4
            # e1 = np.linspace(self.a,quarter_point,half+1)
            # e2 = np.linspace(quarter_point,self.b,N+1-half)[1:]
            # self.x = np.concatenate((e1,e2))
        else:
            xs_cheb = np.array([np.cos(np.pi * j / (self.N)) for j in range(self.N + 1)])
            if np.abs(self.a) == np.abs(self.b) == 1:
                self.x = xs_cheb
            else:
                # range is [a,b]
                self.x = np.array([(self.b - self.a) * x / 2 + (self.a + self.b) / 2 for x in xs_cheb])
        self.h = (b - a) / N
        self.p_prime = p_prime


    #Linear basis functions
    def phi(self, i, x):
        if i != 0 and self.x[i] >= x and x >= self.x[i - 1]:
            return (x - self.x[i - 1]) / (self.x[i]-self.x[i-1])
        elif i != self.N and self.x[i + 1] >= x and x >= self.x[i]:
            return (self.x[i + 1] - x) / (self.x[i+1]-self.x[i])
        else:
            return 0

    def phi_prime(self, i, x):
        if i!= 0 and self.x[i] >= x and x >= self.x[i - 1]:
            return 1 / (self.x[i]-self.x[i-1])
        elif i!= self.N and self.x[i + 1] >= x and x >= self.x[i]:
            return -1 / (self.x[i+1]-self.x[i])
        else:
            return 0

    def fem(self):
        """Solve the SL equation using Finite Element Method
        Cite : An Introduction to Numerical Analysis (Endre Süli, David F. Mayers)
        --------------------------------------------------------------------------
        Returns:
            Array,Array : solution y values, xs
        """
        A = lil_matrix((self.N + 1, self.N + 1))  # instead of lil_matrix
        v = np.zeros(self.N + 1)
        u = np.zeros(self.N + 1)
        for i in range(self.N):
            diag_int1 = lambda x: self.r(x) * self.phi(i, x) ** 2
            diag_int2 = lambda x: self.p(x) * self.phi_prime(i, x) ** 2
            if i == 0:
                diag = (
                    sp.integrate.quad(diag_int1, self.x[i], self.x[i + 1])[0]
                    + sp.integrate.quad(diag_int2, self.x[i], self.x[i + 1])[0]
                )
            elif i == self.N:
                diag = (
                    sp.integrate.quad(diag_int1, self.x[i - 1], self.x[i])[0]
                    + sp.integrate.quad(diag_int2, self.x[i -1], self.x[i])[0]
                )
            else:
                diag = (
                    sp.integrate.quad(diag_int1, self.x[i - 1], self.x[i + 1])[0]
                    + sp.integrate.quad(diag_int2, self.x[i -1], self.x[i + 1])[0]
                )
            
            A[i, i] = diag
            if i > 0:
                lower_int1 = lambda x: self.r(x) * self.phi(i, x) * self.phi(i - 1, x)
                lower_int2 = (
                    lambda x: self.p(x)
                    * self.phi_prime(i, x)
                    * self.phi_prime(i - 1, x)
                )
                lower = (
                    sp.integrate.quad(lower_int1, self.x[i - 1], self.x[i])[0]
                    + sp.integrate.quad(lower_int2, self.x[i - 1], self.x[i])[0]
                )
                A[i, i - 1] = lower
            if i < self.N - 1:
                upper_int1 = lambda x: self.r(x) * self.phi(i, x) * self.phi(i + 1, x)
                upper_int2 = (
                    lambda x: self.p(x)
                    * self.phi_prime(i, x)
                    * self.phi_prime(i + 1, x)
                )
                upper = (
                    sp.integrate.quad(upper_int1, self.x[i], self.x[i + 1])[0]
                    + sp.integrate.quad(upper_int2, self.x[i], self.x[i + 1])[0]
                )
                A[i, i + 1] = upper

            v_int = lambda x: self.phi(i, x) * self.f(x)
            # * bi = <f,phi_i> - <psi,phi_i>
            # * but psi = A * phi_0 + B*phi_n where u(a) = A and u(b) = B so A = B = 0
            v[i] = sp.integrate.quad(v_int, self.x[i - 1], self.x[i + 1])[
                0
            ]  # self.h * self.f(self.x[i])
        A = A.tocsr()
        A = A[1:-1, 1:-1]
        v = v[1:-1]
        u[1:-1] = spsolve(A, v)
        return u, self.x


    #Master element function
    def phi_stand(self, i, xi):
        if i == 0:
            return  1 - xi
        elif i == 1:
            return xi
        else:
            return 0

    def phi_stand_prime(self, i, xi):
        if i == 0:
            return -1
        elif i == 1:
            return 1
        else:
            return 0 

    def fem_stand(self):
        """Solves SL equation using the Master Element method with linear elements 
        Cite :  T. Wick. Numerical methods for partial differential equations. Hannover: Institutionelles
Repositorium der Leibniz Universität Hannover, 2022.
--------------------------------------------------------
Returns: Array,Array : solution y values, xs
"""
        A = lil_matrix((self.N + 1, self.N + 1))  
        v = np.zeros(self.N + 1)
        u = np.zeros(self.N + 1)
        #master element a and b
        ma = 0
        mb = 1
        for i in range(self.N):
            h_i = self.x[i+1]-self.x[i]
            #ξ -> x
            x_transform = lambda xi: h_i * xi + self.x[i]
            x_inv_trans = lambda x: (x-self.x[i])/h_i
            dT = h_i
            dT_inv = 1/h_i
            #jacobian
            J = h_i
            
            for a in [0,1]:
                g_i = i + a
                diag_int1 = (
                    lambda xi: self.r(x_transform(xi)) * self.phi_stand(a, xi) ** 2  * J
                )
                diag_int2 = (
                    lambda xi: self.p(x_transform(xi))
                    * self.phi_stand_prime(a, xi) ** 2 * dT_inv**2 * J 
                )
                A[g_i,g_i] += (sp.integrate.quad(diag_int1,ma,mb)[0] + sp.integrate.quad(diag_int2,ma,mb)[0])
            
            off_int1 = lambda xi: self.r(x_transform(xi)) * self.phi_stand(0, xi) * self.phi_stand(1, xi) * J
            off_int2 = lambda xi: self.p(x_transform(xi)) * self.phi_stand_prime(0, xi) * self.phi_stand_prime(1, xi) * dT_inv**2 * J
            val = sp.integrate.quad(off_int1, ma, mb)[0] + sp.integrate.quad(off_int2, ma, mb)[0]
            A[i, i+1] += val
            A[i+1, i] += val  
            #Instead of the following by symmetry
            # if i > 0:
            #     lower_int1 = (
            #         lambda xi: self.r(x_transform(xi))
            #         * self.phi_stand(1, xi)
            #         * self.phi_stand(0, xi)
            #         * dT_inv**2*J
            #     )
            #     lower_int2 = (
            #         lambda xi: self.p(x_transform(xi))
            #         * self.phi_stand_prime(1, xi)
            #         * self.phi_stand_prime(0, xi)
            #         * dT_inv**2*J
            #     )
            #     lower = (
            #         sp.integrate.quad(lower_int1, ma, mb)[0]
            #         + sp.integrate.quad(lower_int2, ma, mb)[0]
            #     )
            #     A[i, i - 1] += lower
            # if i < self.N - 1:
            #     upper_int1 = (
            #         lambda xi: self.r(x_transform(xi))
            #         * self.phi_stand(0, xi)
            #         * self.phi_stand(1, xi)
            #         * dT_inv**2*J
            #     )
            #     upper_int2 = (
            #         lambda xi: self.p(x_transform(xi))
            #         * self.phi_stand_prime(0, xi)
            #         * self.phi_stand_prime(1, xi)
            #         * dT_inv**2 * J
            #     )
            #     upper = (
            #         sp.integrate.quad(upper_int1, ma, mb)[0]
            #         + sp.integrate.quad(upper_int2, ma, mb)[0]
            #     )
            #     A[i, i + 1] += upper

            v_int = lambda xi: self.phi_stand(0, xi) * self.f(x_transform(xi)) * J
            v_int2 = lambda xi: self.phi_stand(1, xi) * self.f(x_transform(xi)) * J
            v[i] += sp.integrate.quad(v_int, ma, mb)[0]
            v[i+1] += sp.integrate.quad(v_int2, ma, mb)[0]
        A = A.tocsr()
        A = A[1:-1, 1:-1]
        v = v[1:-1]
        u[1:-1] = spsolve(A, v)
        return u, self.x

    def phi_stand_quad(self, i, xi):
        if i == 0:
            return 2*xi**2 - 3*xi + 1
        elif i == 1:
            return 4*xi*(1 - xi)
        elif i == 2:
            return 2*xi**2 - xi
        else:
            return 0

    def phi_stand_quad_prime(self, i, xi):
        if i == 0:
            return 4*xi - 3
        elif i == 1:
            return 4 - 8*xi
        elif i == 2:
            return 4*xi - 1
        else:
            return 0
        
    
    def generate_p2_nodes(self,x_og):
        x_2 = np.zeros(2 * self.N + 1)
        
        # Set original nodes at even indices
        x_2[0::2] = x_og
        # Set midpoints at odd indices
        x_2[1::2] = (x_og[:-1] + x_og[1:]) / 2
        
        return x_2

    def fem_stand_quad(self):
        """Solves SL equation using the Master Element method with quadratic elements elements 
        Cite :  T. Wick. Numerical methods for partial differential equations. Hannover: Institutionelles
        Repositorium der Leibniz Universität Hannover, 2022.
        --------------------------------------------------------
        Returns: Array,Array : solution y values, xs"""
        A = lil_matrix((2*self.N + 1, 2*self.N + 1))
        v = np.zeros(2*self.N + 1)
        u = np.zeros(2*self.N + 1)
        x_og = self.x
        self.x = self.generate_p2_nodes(x_og)

        #master element a and b
        ma = 0
        mb = 1
        for i in range(self.N):
            h_i = self.x[2*i+2] - self.x[2*i]
            x_transform = lambda xi: h_i * xi + self.x[2*i]
            dT_inv = 1/h_i
            J = h_i
            #global elements
            g_i = [2*i, 2*i + 1, 2*i + 2]

            for a in [0, 1, 2]:
                for b in [0, 1, 2]:
                    int1 = lambda xi: self.r(x_transform(xi)) * self.phi_stand_quad(a, xi) * self.phi_stand_quad(b, xi) * J

                    int2 = lambda xi: self.p(x_transform(xi)) * self.phi_stand_quad_prime(a, xi) * self.phi_stand_quad_prime(b, xi) * dT_inv**2 * J

                    val = (
                        sp.integrate.quad(int1, ma, mb)[0] +
                        sp.integrate.quad(int2, ma, mb)[0]
                    )
                    
                    A[g_i[a], g_i[b]] += val

                # RHS
                b_int = lambda xi: self.f(x_transform(xi)) * self.phi_stand_quad(a, xi) * J
                v[g_i[a]] += sp.integrate.quad(b_int, ma, mb)[0]
        A = A.tocsr()
        A = A[1:-1, 1:-1]
        v = v[1:-1]
        u[1:-1] = spsolve(A, v)
        return u, self.x
    
    #Sth that I didn't finish implementing but should be kinda correct
    def phi_quart(self, i, x):
        # xl and xr are the midpoints to the right and left
        if i == 0:
            if x <= self.x[i + 1]:
                xr = (self.x[i + 1] + self.x[i]) / 2
                return (
                    (x - self.x[i + 1])
                    * (x - xr)
                    / ((self.x[i] - xr) * (self.x[i] - self.x[i + 1]))
                )
            else:
                return 0

        elif i == self.N:
            if self.x[i - 1] <= x:
                xl = (self.x[i - 1] + self.x[i]) / 2
                return (
                    (x - self.x[i - 1])
                    * (x - xl)
                    / ((self.x[i] - xl) * (self.x[i] - self.x[i - 1]))
                )
            else:
                return 0
        else:
            if self.x[i - 1] <= x <= self.x[i + 1]:
                xl = (self.x[i - 1] + self.x[i]) / 2
                xr = (self.x[i + 1] + self.x[i]) / 2

                return (
                    (x - self.x[i - 1])
                    * (x - xl)
                    * (x - xr)
                    * (x - self.x[i + 1])
                    / (
                        (self.x[i] - xl)
                        * (self.x[i] - xr)
                        * (self.x[i] - self.x[i - 1])
                        * (self.x[i] - self.x[i + 1])
                    )
                )
            else:
                return 0

    def phi_quart_prime(self, i, x):
        # xl and xr are the midpoints to the right and left
        if i == 0:
            if x <= self.x[i + 1]:
                xr = (self.x[i + 1] + self.x[i]) / 2
                return ((x - self.x[i + 1]) + (x - xr)) / (
                    (self.x[i] - xr) * (self.x[i] - self.x[i + 1])
                )
            else:
                return 0

        elif i == self.N:
            if self.x[i - 1] <= x:
                xl = (self.x[i - 1] + self.x[i]) / 2
                return ((x - self.x[i - 1]) + (x - xl)) / (
                    (self.x[i] - xl) * (self.x[i] - self.x[i - 1])
                )
            else:
                return 0
        else:
            if self.x[i - 1] <= x <= self.x[i + 1]:
                xl = (self.x[i - 1] + self.x[i]) / 2
                xr = (self.x[i + 1] + self.x[i]) / 2

                a = self.x[i - 1]
                b = xl
                c = xr
                d = self.x[i + 1]
                num = (
                    (x - b) * (x - c) * (x - d)
                    + (x - a) * (x - c) * (x - d)
                    + (x - a) * (x - b) * (x - d)
                    + (x - a) * (x - b) * (x - c)
                )
                den = (
                    (self.x[i] - a)
                    * (self.x[i] - b)
                    * (self.x[i] - c)
                    * (self.x[i] - d)
                )

                return num / den
            else:
                return 0

    def fem_quart(self):
        """Solve the SL equation using Finite Element Method (doesn't really work tho)
        Cite : An Introduction to Numerical Analysis (Endre Süli, David F. Mayers)
        Returns:
            Array: solution y values
        """
        A = csr_matrix((self.N + 1, self.N + 1))  # instead of lil_matrix
        v = np.zeros(self.N + 1)
        u = np.zeros(self.N + 1)
        for i in range(self.N):
            diag_int1 = lambda x: self.r(x) * self.phi_quart(i, x) ** 2
            diag_int2 = lambda x: self.p(x) * self.phi_quart_prime(i, x) ** 2

            diag = (
                sp.integrate.quad(diag_int1, self.x[i - 1], self.x[i + 1])[0]
                + sp.integrate.quad(diag_int2, self.x[i - 1], self.x[i + 1])[0]
            )
            A[i, i] = diag
            if i > 0:
                lower_int1 = (
                    lambda x: self.r(x)
                    * self.phi_quart(i, x)
                    * self.phi_quart(i - 1, x)
                )
                lower_int2 = (
                    lambda x: self.p(x)
                    * self.phi_prime(i, x)
                    * self.phi_prime(i - 1, x)
                )
                lower = (
                    sp.integrate.quad(lower_int1, self.x[i - 1], self.x[i])[0]
                    + sp.integrate.quad(lower_int2, self.x[i - 1], self.x[i])[0]
                )
                A[i, i - 1] = lower
            if i < self.N - 1:
                upper_int1 = (
                    lambda x: self.r(x)
                    * self.phi_quart(i, x)
                    * self.phi_quart(i + 1, x)
                )
                upper_int2 = (
                    lambda x: self.p(x)
                    * self.phi_prime(i, x)
                    * self.phi_prime(i + 1, x)
                )
                upper = (
                    sp.integrate.quad(upper_int1, self.x[i], self.x[i + 1])[0]
                    + sp.integrate.quad(upper_int2, self.x[i], self.x[i + 1])[0]
                )
                A[i, i + 1] = upper

            v_int = lambda x: self.phi_quart(i, x) * self.f(x)
            # * bi = <f,phi_i> - <psi,phi_i>
            # * but psi = A * phi_0 + B*phi_n where u(a) = A and u(b) = B so A = B = 0
            v[i] = sp.integrate.quad(v_int, self.x[i - 1], self.x[i + 1])[
                0
            ]  # self.h * self.f(self.x[i])

        A = A[1:-1, 1:-1]
        v = v[1:-1]
        u[1:-1] = spsolve(A, v)
        return u, self.x

    def chebyshev_diff_matrix(self, xs):
        """Chebyshev differentiation matrix:
        Cite: Trefethen, Lloyd N. — Spectral Methods in MATLAB
        """
        # as shown in book p.54
        c = np.ones(self.N + 1)
        c[0] = 2
        c[-1] = 2
        c = c * (-1) ** np.arange(self.N + 1)  # c = [2,1,1,...,1,2]

        X = np.tile(xs, (self.N + 1, 1))  # [xs,xs,..,xs]
        dX = X - X.T  # diag 0
        D = (np.outer(c, 1 / c)) / (dX + np.eye(self.N + 1))  # eye diag 1s else 0
        D = D - np.diag(np.sum(D, axis=1))
        D = -1 * D
        return D

    def sm(self):
        """Solve the SL equation using Spectral Element Method
        Cite:  Trefethen, Lloyd N. — Spectral Methods in MATLAB
        Returns:
           list,list: solution, x values
        """

        xs_cheb = np.array([np.cos(np.pi * j / (self.N)) for j in range(self.N + 1)])
        # Collocation points - Gauss-Lobatto
        if np.abs(self.a) == np.abs(self.b) == 1:
            xs = xs_cheb
        else:
            # scales cheb points to [a,b]
            xs = np.array([(self.b - self.a) * x / 2 + (self.a + self.b) / 2 for x in xs_cheb])

        D = self.chebyshev_diff_matrix(xs)
        D2 = D @ D

        # Au = f
        f = np.array([self.f(x) for x in xs])
        p_vals = np.diag([self.p(x) for x in xs])
        r_vals = np.diag([self.r(x) for x in xs])
        if self.p_prime == None:
            # calculate derivative using Cheb matrix
            p_prime = D @ np.array([self.p(x) for x in xs])
        else:
            p_prime = np.array([self.p_prime(x) for x in xs])

        p_prime_vals = np.diag(p_prime)

        A = -(p_vals @ D2 + p_prime_vals @ D) + r_vals

        if np.allclose(f, 0):
            print(
                f"WARNING: This will most likely return 0, recommend solving the Eigenproblem instead"
            )
        u = np.zeros(self.N + 1)
        A = A[1:-1, 1:-1]
        f = f[1:-1]
        u[1:-1] = np.linalg.solve(A, f)
        return u, xs

    def differentiate(self, f):
        xs = np.array([np.cos(np.pi * j / (self.N)) for j in range(self.N + 1)])
        f_vals = np.array([f(x) for x in xs])
        D = self.chebyshev_diff_matrix(xs)

        f_prime = D @ f_vals
        return f_prime, xs
    
    

    def sm2(self):
        N = self.N
        xs = self.legendre_gll_nodes(N)
        w = self.gll_weights(N,xs)

        # Derivative matrix
        D = self.deriv_matrix_local(N,xs) 
        # Evaluate functions
        p_vals = np.array([self.p(x) for x in xs])
        r_vals = np.array([self.r(x) for x in xs])
        f_vals = np.array([self.f(x) for x in xs])

        # Stiffness + Mass matrix
        
        
        A = np.diag(w * r_vals) + (D.T @ np.diag(w * p_vals) @ D)
        #Can do ^ instead of below 
        #A = lil_matrix((N + 1, N + 1))
        # for i in range(N + 1):
        #     for j in range(N + 1):
        #         A[i, j] = np.sum(w * p_vals * D[:, j] * D[:, i]) + w[i] * r_vals[i] * (i == j)

        v = w * f_vals

        #Boundary condition
        
        A_inner = A[1:-1, 1:-1]
        v_inner = v[1:-1]
        u = np.zeros(N + 1)
        u[1:-1] = spsolve(A_inner, v_inner)

        return u, xs
    
    def legendre_gll_nodes(self,N):
        #Spectral/hp Elem methods for CFD Karniadakis
        if N == 1:
            return np.array([-1.0, 1.0])

        # Use Chebyshev-Gauss-Lobatto as initial guess for Newton like Karniadakis
        x = np.cos(np.pi * np.arange(N + 1) / N)

        # Newton iteration to find roots of (1 - x^2) * L_N'(x) = 0
        L = sp.special.legendre(N)
        Lp = L.deriv()
        
        for _ in range(100):
            dx = -((1 - x**2) * Lp(x)) / (-2 * x * Lp(x) + (1 - x**2) * Lp.deriv()(x))
            x += dx
            if np.max(np.abs(dx)) < 1e-14:
                break

        return np.sort(x)
    
    def gll_weights(self, N,xs):
        #Boyd 
        L = sp.special.legendre(N)
        weights = 2 / (N * (N + 1)) / (L(xs)**2)
        weights[0] = 2 / (N * (N + 1))  # Endpoint weights
        weights[-1] = 2 / (N * (N + 1))
        return weights
    
    def deriv_matrix_local(self,N,xj):
        #Boyd
        L_n = sp.special.legendre(N)
        L_n_prime = L_n.deriv()
        
        D = np.zeros((N+1,N+1))
        
        #note those 2 can be different like x_evaluate points and x_collocation nodes
        for i in range(N+1):
            for p in range(N+1):
                nom = (1-xj[i])**2*L_n_prime(xj[i])
                den = N*(N+1) * L_n(xj[p])
                if i == p == 0:
                    D[i,p] = -N*(N+1)/4
                elif i == p == N:
                    D[i,p] = N*(N+1)/4
                 
                elif i != p:       
                    D[i,p] = L_n(xj[i])/(L_n(xj[p])*(xj[i]-xj[p]))#nom/ (den*(xj[i]-xj[p]))
                else:
                    D[i,p] = 0
        
        return D
    
    
    def sem(self, ne,N=None):
        '''A bit from Citation for published version (APA):
Timmermans, L. J. P., Jansen, J. K. M., & Vosse, van de, F. N. (1990).
A description of the fundamentals of the
spectral element method. (DCT rapporten; Vol. 1990.041). Technische Universiteit Eindhoven'''
        #so it's not the same as #num elements for FEM and SEM
        #Here N = Polynomial degree
        if N == None:
            N = self.N

        #Get Gauss-Lobbato-Legendre nodes degree N
        xs = self.legendre_gll_nodes(N)
        ws = self.gll_weights(N,xs)

        #ne number of elements between a and b which are -1 and 1 in this case
        element_nodes = np.linspace(self.a, self.b, ne+1)
        #Uncomment below if you want half the elements to be in the first quater of the interval
        # half = ne // 2
        # quarter_point = (3*self.a+self.b)/4
        # e1 = np.linspace(self.a,quarter_point,half+1)
        # e2 = np.linspace(quarter_point,self.b,ne+1-half)[1:]
        # element_nodes = np.concatenate((e1,e2))
        
        l = np.diff(element_nodes) #len of each interval same as before x[i+1]-x[i] just a vector now
        ng = ne * N + 1  # (L-1)*N+N+1 -igel # of global nodes
        
        #Initialize global matrices
        A_global = lil_matrix((ng, ng))
        b_global = np.zeros(ng)
        u = np.zeros(ng)
        x_g = np.zeros(ng) #keep track of global integration nodes
        for j in range(ne):
            J = l[j]/2
            Ji = 2/l[j]
            
            #local nodes and weights on element j
            xj = [element_nodes[j] +  (xs[i]+1)*J for i in range(N+1)]
            wj = ws*J
            
            #Keep track of global nodes
            #Old way need x_g = []
            # if j == 0:
            #     x_g.extend(xj)
            # else:
            #     x_g.extend(xj[:-1])
            
            #New way so x_g can be a np.array
            for i in range(N + 1):
                global_index = j * N + i
                x_g[global_index] = xj[i]
            #print(x_g) 
            #Derivative matrix
            #D_loc = self.deriv_matrix_local(N,xj) 
            D_loc = self.deriv_matrix_local(N,xs)  
            D_loc_scaled = D_loc *Ji

            # Function values on element j
            p_vals = np.array([self.p(x) for x in xj])
            r_vals = np.array([self.r(x) for x in xj])
            f_vals = np.array([self.f(x) for x in xj])

            # Build local matrices
            A_loc = (D_loc_scaled.T @ np.diag(wj * p_vals) @ D_loc_scaled)+ np.diag(wj * r_vals)
            b_loc = wj * f_vals
            #Above is a shorter version nof the below
            # A_loc = np.zeros((N + 1, N + 1))
            # b_loc = np.zeros(N + 1)
            # for i in range(N + 1):
            #     for p in range(N + 1):
            #         A_loc[i, p] = np.sum(wj * p_vals * D_loc[:,p] * D_loc[:,i]) + wj[i]*r_vals[i]*(i==p)

            #     b_loc[i] = wj[i] * f_vals[i]
           
            
            # Map local to global and assemble
            g_nodes =  [j*N + i for i in range(N+1)]
            for i in range(N + 1):
                gi = g_nodes[i]
                b_global[gi] += b_loc[i]
                for e in range(N + 1):
                    ge = g_nodes[e]
                    A_global[gi, ge] += A_loc[i, e]
                    
        # Apply Dirichlet BCs at endpoints i.e. matrix condensation
        A_global = A_global.tocsr()
        A_global = A_global[1:-1,1:-1]
        b_global = b_global[1:-1]
        u[1:-1] = spsolve(A_global, b_global)
        #Element midpoints and averages 
        elem_midpoints = []
        u_avg = []
        for i in range(ne):
            elem_midpoints.append((element_nodes[i]+element_nodes[i+1])/2)
            start = i*N
            end = start+N+1
            u_avg.append(np.mean(u[start:end]))
        return u,x_g
