import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix, lil_matrix


class EigenvalueProblem:

    def __init__(self, p, r, a, b, N, p_prime):
        """Euqation -d/dx(p(x)u'(x))+r(x)u(x)=0"""
        self.p = p
        self.r = r

        if a <= b:
            self.a = a
            self.b = b
        else:
            self.a = b
            self.b = a

        self.N = N
        self.h = (b - a) / N
        self.p_prime = p_prime
        self.x = np.linspace(a, b, N + 1)

    def phi_stand(self, i, xi):
        if i == 0:
            return 1 - xi
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
        A = csr_matrix((self.N + 1, self.N + 1))
        B = csr_matrix((self.N + 1, self.N + 1))
        # master element a and b
        ma = 0
        mb = 1
        for i in range(self.N):
            h_i = self.x[i + 1] - self.x[i]
            # ξ -> x
            x_transform = lambda xi: h_i * xi + self.x[i]
            dT_inv = 1 / h_i
            # jacobian
            J = h_i

            for a in [0, 1]:
                g_i = i + a
                b_diag = (
                    lambda xi: self.r(x_transform(xi)) * self.phi_stand(a, xi) ** 2 * J
                )
                a_diag = (
                    lambda xi: self.p(x_transform(xi))
                    * self.phi_stand_prime(a, xi) ** 2
                    * dT_inv**2
                    * J
                )
                A[g_i, g_i] += sp.integrate.quad(a_diag, ma, mb)[0]
                B[g_i, g_i] += sp.integrate.quad(b_diag, ma, mb)[0]

            b_off = (
                lambda xi: self.r(x_transform(xi))
                * self.phi_stand(0, xi)
                * self.phi_stand(1, xi)
                * J
            )
            a_off = (
                lambda xi: self.p(x_transform(xi))
                * self.phi_stand_prime(0, xi)
                * self.phi_stand_prime(1, xi)
                * dT_inv**2
                * J
            )
            a_off_val = sp.integrate.quad(a_off, ma, mb)[0]
            b_off_val = sp.integrate.quad(b_off, ma, mb)[0]
            A[i, i + 1] += a_off_val
            A[i + 1, i] += a_off_val

            B[i, i + 1] += b_off_val
            B[i + 1, i] += b_off_val

        A = A[1:-1, 1:-1]
        B = B[1:-1, 1:-1]
        l, v = sp.linalg.eig(A.toarray(), B.toarray())
        # idx = np.argsort(l)
        # l = l[idx]
        # v = v[:, idx]
        return l, v, self.x

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

        l = np.diff(element_nodes) #len of each interval same as before x[i+1]-x[i] just a vector now
        ng = ne * N + 1  # (L-1)*N+N+1 -igel # of global nodes
        
        #Initialize global matrices
        A_global = lil_matrix((ng, ng))
        B_global = lil_matrix((ng, ng))
        u = np.zeros(ng)
        x_g = np.zeros(ng) #keep track of global integration nodes
        for j in range(ne):
            J = l[j]/2
            Ji = 2/l[j]
            
            #local nodes and weights on element j
            xj = [element_nodes[j] +  (xs[i]+1)*J for i in range(N+1)]           
            wj = ws*J

            
            #New way so x_g can be a np.array
            for i in range(N + 1):
                global_index = j * N + i
                x_g[global_index] = xj[i]
            
            #Derivative matrix
            D_loc = self.deriv_matrix_local(N,xs)  
            D_loc_scaled = D_loc *Ji

            # Function values on element j
            p_vals = np.array([self.p(x) for x in xj])
            r_vals = np.array([self.r(x) for x in xj])

            # Build local matrices
            A_loc = (D_loc_scaled.T @ np.diag(wj * p_vals) @ D_loc_scaled)
            B_loc = np.diag(wj * r_vals)
            
            # Map local to global and assemble
            g_nodes =  [j*N + i for i in range(N+1)]
            for i in range(N + 1):
                gi = g_nodes[i]
                for e in range(N + 1):
                    ge = g_nodes[e]
                    A_global[gi, ge] += A_loc[i, e]
                    B_global[gi, ge] += B_loc[i, e]
                    
        # Apply Dirichlet BCs at endpoints i.e. matrix condensation
        A_global = A_global.tocsr()
        A_global = A_global[1:-1,1:-1]
        B_global = B_global.tocsr()
        B_global = B_global[1:-1,1:-1]
        l,v = sp.linalg.eig(A_global.toarray(),B_global.toarray())
        #Element midpoints and averages 

        return l,v,x_g