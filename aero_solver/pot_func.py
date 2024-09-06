import math
import scipy.integrate as inte
import numpy as np
import numpy as np

PI = math.pi
PI_inv = 1 / math.pi

class aero_solver_osc_flat:
    def __init__(self, kin, U_ref, t_step, alpha_eff):
        self.kin = kin
        self.U_ref = U_ref

        self.v_core = 1.3*t_step*U_ref
        self.alpha_eff = alpha_eff

    def xin2body(self, x, t):
        return x + self.U_ref*t

    def yin2body(self, y, t):
        return y - self.kin.h(t)

    def bodyin2x(self, x,t):
        return x - self.U_ref*t

    def bodyin2y(self, y, t):
        return y + self.kin.h(t)


    def g_trans(self, theta, c):
        return 0.5 * c * (1 - np.cos(theta))

    def dphideta(self, xi_n, eta_n , Gamma_n):

        const = 0.5 * Gamma_n * PI_inv

        func = lambda xi:  const * (
            (-eta_n * math.sin(self.alpha_eff) - (xi - xi_n) * math.cos(self.alpha_eff))  
            / 
            np.sqrt((eta_n**2 + (xi - xi_n)**2)**2 + self.v_core**4)
                             )
        return func

    def W_0(self, t):

        return lambda xi: - self.U_ref*math.sin(self.alpha_eff) + (self.kin.h_dot(t)) * math.cos(self.alpha_eff)

    def W_0_fast_1(self, t):

        return - self.U_ref*math.sin(self.alpha_eff) + (self.kin.h_dot(t)) * math.cos(self.alpha_eff)


    def V_ind_b_fast_2(self, gamma, xi_n, eta_n, c):
        '''
        - takes in vorticity distribution gamma
        - find the induced velocity of the vorticity distribution at point (xi_n, eta_n)
        '''

        x = np.linspace(0.0001, c, 513, endpoint=True)

        integrand_u = lambda xi: gamma(xi) * (0.0 - eta_n) / ((((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)**0.5)

        def_int_u = inte.trapezoid(integrand_u(x),x)

        u_ind = 0.5 * PI_inv * def_int_u

        integrand_v = lambda xi: gamma(xi) *  (xi -  xi_n) / ((((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)**0.5)

        def_int_v = inte.trapezoid(integrand_v(x),x)

        v_ind = -0.5 * PI_inv * def_int_v 

        return u_ind, v_ind

    def V_ind_ub(self, xi, eta, xi_n, eta_n, gamma):

        '''
        - calculates the induced velocity at a point by another vortex blob centred at (xi_n, eta_n)
        '''

        u_ind_ub = 0.5 * gamma * PI_inv * (eta - eta_n) / math.sqrt(((xi - xi_n)**2 + (eta - eta_n)**2)**2 + self.v_core**4)

        v_ind_ub = -0.5 * gamma * PI_inv * (xi - xi_n) / math.sqrt(((xi - xi_n)**2 + (eta - eta_n)**2)**2 + self.v_core**4)

        return u_ind_ub, v_ind_ub
    


    def V_ind_tot_field(self, x1_N, y1_N, x2_N, y2_N, Gamma_N, fourier, no_gamma, U, c, t):
            
        u_ind, v_ind = self.V_ind_ub_field(x1_N, y1_N, x2_N, y2_N, Gamma_N, no_gamma)

        # Finding induced velocity at each vortex
        for n in range(len(u_ind)):
            
            # Induced velocity on a vortex by the bounded vortex sheet            
            trans = lambda xi: np.arccos(1 - 2*xi/c)
            gamma = lambda xi: 2* U * (fourier[0] * (1 + np.cos(trans(xi)))/np.sin(trans(xi)) + fourier[1] * np.sin(trans(xi))) + fourier[2] * np.sin(2*trans(xi)) + fourier[3] * np.sin(3*trans(xi)) #+ fourier[4] * np.sin(4*trans(xi)) + fourier[5] * np.sin(5*trans(xi))

            u_ind_p, v_ind_p = self.V_ind_b_fast_2(gamma, self.xin2body(x1_N[n],t), self.yin2body(y1_N[n],t), c)

            u_ind[n] += u_ind_p #+ U_ref
            v_ind[n] += v_ind_p #+ pot.hdot(t) 

        
        return u_ind, v_ind
    
    
    
    def V_ind_ub_field(self, x1_N, y1_N, x2_N, y2_N, Gamma_N,dd):
        '''
        calculates induced velocity at (x1,y1) by vortices at (x2,y2)
        '''
        # Reshape the input arrays to enable broadcasting
        x1_N = np.reshape(x1_N, (-1, 1))  # Shape: (len(x1_N), 1)
        y1_N = np.reshape(y1_N, (-1, 1))  # Shape: (len(y1_N), 1)
        # Gamma_N = np.reshape(Gamma_N, (-1, 1))  # Shape: (len(y1_N), 1)

        # Compute the difference arrays (x1_N - x2_N) and (y1_N - y2_N)
        dx = x1_N - x2_N  # Shape: (len(x1_N), len(x2_N))
        dy = y1_N - y2_N  # Shape: (len(y1_N), len(y2_N))

        # Calculate the squared distance and avoid division by zero by adding a small value (epsilon)
        r_squared = dx**2 + dy**2
        epsilon = 1e-10  # Small value to prevent division by zero
        r_squared = np.where(r_squared == 0, epsilon, r_squared)

        # Calculate induced velocities using broadcasting
        u_ind_p =  dy*( 1 / (2 * np.pi * (r_squared**2 + self.v_core**4)**0.5)) * Gamma_N # Shape: (len(x1_N), len(x2_N))
        v_ind_p =  dx*(-1 / (2 * np.pi * (r_squared**2 + self.v_core**4)**0.5)) * Gamma_N # Shape: (len(y1_N), len(y2_N))

        # Sum along the second axis to accumulate the induced velocities
        u_ind = np.sum(u_ind_p, axis=1)
        v_ind = np.sum(v_ind_p, axis=1)

        return u_ind, v_ind


    def fourier_gamma_calc(self, A_no, Gamma_N, eta_N, xi_N, N, c, t):

        x = np.linspace(0.0, np.pi, 513, endpoint=True)

        fourier = np.zeros(A_no)
        U_ref_inv = 1 / self.U_ref

        g_trans = lambda theta: 0.5*c*(1-np.cos(theta))

        # Computing Fourier coefficients
        for i in range(len(fourier)):
            fourier[i] = 0.0
            # Computing A_0 in fourier series of vorticity distribution on the bound vortex
            if i == 0: 
                if N == 0: # solving for t = 0
                    fourier[i] = - PI_inv * U_ref_inv * self.W_0_fast_1(t) * np.pi
                else: # solving for t > 0

                    fourier[i] = self.W_0_fast_1(t) * np.pi

                    for n in range(N):                            
                        Gamma_n = Gamma_N[n]
                        eta_n   = eta_N[n]
                        xi_n    = xi_N[n]
                        dphideta = self.dphideta(xi_n, eta_n, Gamma_n)
                        integrand_n = lambda theta: dphideta(g_trans(theta))

                        A_int = inte.trapezoid(integrand_n(x), x)

                        fourier[i] -= A_int
                    fourier[i] *= - 1.0 / np.pi / self.U_ref
            # Computing A_n in fourier series of vorticity distribution on the bound vortex
            else:
                if N == 0: # solving for t = 0
                    integrand_n = lambda theta: self.W_0_fast_1(t) * math.cos(i * theta)
                    fourier[i], extra = inte.quad(integrand_n , 0.0, np.pi)
                    fourier[i] *= 2.0 / np.pi / self.U_ref
                else: # solving for t > 0
                    integrand_n = lambda theta: self.W_0_fast_1(t) * np.cos(i * theta)
                    fourier[i], extra = inte.quad(integrand_n , 0.0, np.pi)
                    for n in range(N):

                        Gamma_n = Gamma_N[n]
                        eta_n   = eta_N[n]
                        xi_n    = xi_N[n]
                        dphideta = self.dphideta(xi_n, eta_n, Gamma_n)
                        integrand_n = lambda theta: dphideta(g_trans(theta)) * np.cos(i * theta)

                        A_int = inte.trapezoid(integrand_n(x), x)

                        fourier[i] -= A_int

                    fourier[i] *= 2.0 / np.pi / self.U_ref

        Gamma_b = np.pi * c * self.U_ref * (fourier[0] + fourier[1] * 0.5)

        return fourier, sum(Gamma_N), Gamma_b + sum(Gamma_N)
    

    # def V_ind_b(self, gamma, xi_n, eta_n):
    #     '''
    #     - takes in vorticity distribution gamma
    #     - find the induced velocity of the vorticity distribution at point (xi_n, eta_n)
    #     '''

    #     integrand_u = lambda xi: gamma(xi) * (0.0 - eta_n) / ((((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)**0.5)

    #     def_int_u, extra = inte.quad(integrand_u, 0, self.c)

    #     u_ind = 0.5 * PI_inv * def_int_u

    #     integrand_v = lambda xi: gamma(xi) * (xi - xi_n) / ((((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)**0.5)

    #     def_int_v, extra = inte.quad(integrand_v, 0, self.c)

    #     v_ind = -0.5 * PI_inv * def_int_v 

    #     return u_ind, v_ind  