import math
import scipy.integrate as inte
import numpy as np
import numpy as np

PI = math.pi
PI_inv = 1 / math.pi


# def hdot(t):
#     if t < 0.5:
#         return 0
#     elif t < 2.5:
#         return 0.5*np.pi*np.sin(0.5*np.pi*(t-0.5))
#     elif t < 3.5:
#         return - np.pi*2*np.sin(np.pi*(t-2.5))
#     elif t < 5.5:
#         return 0.5*np.pi*np.sin(0.5*np.pi*(t-3.5))
    
#     elif t < 6:
#         return 0
#     elif t < 8:
#         return 0.5*np.pi*np.sin(0.5*np.pi*(t-6))
#     elif t < 9:
#         return - np.pi*2*np.sin(np.pi*(t-8))
#     elif t < 11:
#         return 0.5*np.pi*np.sin(0.5*np.pi*(t-9))

# def h(t):
#     if t < 0.5:
#         return 0
#     elif t < 2.5:
#         return 1 - np.cos(0.5*np.pi*(t-0.5))
#     elif t < 3.5:
#         return 2*np.cos(np.pi*(t-2.5))
#     elif t < 5.5:
#         return - 1 - np.cos(0.5*np.pi*(t-3.5))
    
#     elif t < 6:
#         return 0
#     elif t < 8:
#         return 1 - np.cos(0.5*np.pi*(t-6))
#     elif t < 9:
#         return 2*np.cos(np.pi*(t-8))
#     elif t < 11:
#         return - 1 - np.cos(0.5*np.pi*(t-9))


class aero_solver_osc_flat:
    def __init__(self, h, hdot, c, U_ref, t_step, alpha_eff):
        self.c = c
        self.h = h
        self.hdot = hdot
        self.U_ref = U_ref

        self.v_core = 1.3*t_step*U_ref
        self.alpha_eff = alpha_eff

    def xin2body(self, x, t):
        return x + self.U_ref*t

    def yin2body(self, y, t):
        return y - self.h(t)

    def bodyin2x(self, x,t, U_ref):
        return x - U_ref*t

    def bodyin2y(self, y, t):
        return y + self.h(t)


    def g_trans(self, theta):
        return 0.5 * self.c * (1 - np.cos(theta))

    def dphideta(self, xi_n, eta_n , Gamma_n):

        const = 0.5 * Gamma_n * PI_inv

        eta = lambda xi: 0.0

        func = lambda xi:  const * (
            (-eta_n * math.sin(self.alpha_eff) - (xi - xi_n) * math.cos(self.alpha_eff))  
            / 
            np.sqrt((eta_n**2 + (xi - xi_n)**2)**2 + self.v_core**4)
                             )
        return func

    def W_0(self, t):

        A = 1
        p = 0.5

        eta = lambda xi: 0.0

        return lambda xi: - self.U_ref*math.sin(self.alpha_eff) + (self.hdot(t)) * math.cos(self.alpha_eff)

    def W_0_fast_1(self, t):

        return - self.U_ref*math.sin(self.alpha_eff) + (self.hdot(t)) * math.cos(self.alpha_eff)

    def V_ind_b(self, gamma, xi_n, eta_n):
        '''
        - takes in vorticity distribution gamma
        - find the induced velocity of the vorticity distribution at point (xi_n, eta_n)
        '''

        integrand_u = lambda xi: gamma(xi) * (0.0 - eta_n) / ((((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)**0.5)

        def_int_u, extra = inte.quad(integrand_u, 0, self.c)

        u_ind = 0.5 * PI_inv * def_int_u

        integrand_v = lambda xi: gamma(xi) * (xi - xi_n) / ((((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)**0.5)

        def_int_v, extra = inte.quad(integrand_v, 0, self.c)

        v_ind = -0.5 * PI_inv * def_int_v 

        return u_ind, v_ind  

    def V_ind_b_fast_2(self, gamma, xi_n, eta_n):
        '''
        - takes in vorticity distribution gamma
        - find the induced velocity of the vorticity distribution at point (xi_n, eta_n)
        '''

        x = np.linspace(0.0001, self.c, 513, endpoint=True)

        integrand_u = lambda xi: gamma(xi) * (0.0 - eta_n) / ((((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)**0.5)

        def_int_u = inte.trapz(integrand_u(x),x)

        u_ind = 0.5 * PI_inv * def_int_u

        integrand_v = lambda xi: gamma(xi) * (xi - xi_n) / ((((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)**0.5)

        def_int_v = inte.trapz(integrand_v(x),x)

        v_ind = -0.5 * PI_inv * def_int_v 

        return u_ind, v_ind

    def V_ind_ub(self, xi, eta, xi_n, eta_n, gamma):

        '''
        - calculates the induced velocity at a point by another vortex blob
        '''

        u_ind_ub = 0.5 * gamma * PI_inv * (eta - eta_n) / math.sqrt(((xi - xi_n)**2 + (eta - eta_n)**2)**2 + self.v_core**4)

        v_ind_ub = -0.5 * gamma * PI_inv * (xi - xi_n) / math.sqrt(((xi - xi_n)**2 + (eta - eta_n)**2)**2 + self.v_core**4)

        return u_ind_ub, v_ind_ub

    def fourier_gamma_calc(self, A_no, Gamma_N, eta_N, xi_N, N, t):

        x = np.linspace(0.0, np.pi, 513, endpoint=True)

        A = np.zeros(A_no)
        U_ref_inv = 1 / self.U_ref

        # Computing Fourier coefficients
        for i in range(len(A)):
            A[i] = 0.0
            # Computing A_0 in fourier series of vorticity distribution on the bound vortex
            if i == 0: 
                if N == 0: # solving for t = 0
                    A[i] = - PI_inv * U_ref_inv * self.W_0_fast_1(t) * self.c
                else: # solving for t > 0

                    A[i] = self.W_0_fast_1(t) * self.c

                    for n in range(N):                            
                        Gamma_n = Gamma_N[n]
                        eta_n   = eta_N[n]
                        xi_n    = xi_N[n]
                        dphideta = self.dphideta(xi_n, eta_n, Gamma_n)
                        integrand_n = lambda theta: dphideta(self.g_trans(theta))

                        A_int = inte.trapz(integrand_n(x), x)

                        A[i] -= A_int
                    A[i] *= - 1.0 / np.pi / self.U_ref
            # Computing A_n in fourier series of vorticity distribution on the bound vortex
            else:
                if N == 0: # solving for t = 0
                    integrand_n = lambda theta: self.W_0_fast_1(t) * math.cos(i * theta)
                    A[i], extra = inte.quad(integrand_n , 0.0, np.pi)
                    A[i] *= 2.0 / np.pi / self.U_ref
                else: # solving for t > 0
                    integrand_n = lambda theta: self.W_0_fast_1(t) * np.cos(i * theta)
                    A[i], extra = inte.quad(integrand_n , 0.0, np.pi)
                    for n in range(N):

                        Gamma_n = Gamma_N[n]
                        eta_n   = eta_N[n]
                        xi_n    = xi_N[n]
                        dphideta = self.dphideta(xi_n, eta_n, Gamma_n)
                        integrand_n = lambda theta: dphideta(self.g_trans(theta)) * np.cos(i * theta)

                        A_int = inte.trapz(integrand_n(x), x)

                        A[i] -= A_int

                    A[i] *= 2.0 / np.pi / self.U_ref

        Gamma_b = np.pi * self.c * self.U_ref * (A[0] + A[1] * 0.5)

        return A, sum(Gamma_N), Gamma_b + sum(Gamma_N)
