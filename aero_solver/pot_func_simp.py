import math
import scipy.integrate as inte
import numpy as np
import numpy as np

PI = math.pi
PI_inv = 1 / math.pi

class aero_solver_osc_flat:
    def __init__(self, kin, U_ref, t_step, alpha_eff,chord):
        self.kin = kin
        self.U_ref_M = U_ref
        self.U_ref = U_ref*math.cos(alpha_eff)
        self.V_ref = U_ref*math.sin(alpha_eff)

        self.v_core = 0.02*chord
        self.alpha_eff = alpha_eff

    def xin2body(self, x, t):
        # return x + self.U_ref*t

        return x + self.kin.pos(t)

    def yin2body(self, y, t):
        return y - self.kin.h(t)

    def bodyin2x(self, x,t):
        # return x - self.U_ref*t
        return x - self.kin.pos(t)

    def bodyin2y(self, y, t):
        return y + self.kin.h(t)

 

    def g_trans(self, theta, c):
        return 0.5 * c * (1 - np.cos(theta))

    def dphideta(self, xi_n, eta_n , Gamma_n):

        const = 0.5 * Gamma_n * PI_inv

        func = lambda xi:  const * (
            (-(xi - xi_n))  
            / 
            np.sqrt((eta_n**2 + (xi - xi_n)**2)**2 + self.v_core**4)
                             )
        return func

    def W_0(self, t):

        return lambda xi:  (self.kin.h_dot(t))

    def W_0_fast_1(self, t):

        return (self.kin.h_dot(t))

    
    def V_ind_b_fast_4(self, fourier, xi_n, eta_n, c, t):

        theta = np.linspace(0,np.pi,100,endpoint=True)

        xi = 0.5 * c * (1 - np.cos(theta))
    
        fourier_inf = fourier[0]*(1+np.cos(theta))

        for i in range(1,len(fourier)):

            fourier_inf += fourier[i]*np.sin(i*theta)*np.sin(theta)

        xis = (xi_n - xi)/np.sqrt(((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)

        etas = (eta_n)   /np.sqrt(((xi_n - xi)**2 + eta_n**2)**2 + self.v_core**4)

        
        u_ind = 0.5 * c * self.U_ref * inte.trapezoid(etas*fourier_inf,theta) / np.pi

        v_ind =-0.5 * c * self.U_ref * inte.trapezoid(xis*fourier_inf,theta) / np.pi

        return u_ind, v_ind

    
    def V_ind_tot_field(self, x1_N, y1_N, x2_N, y2_N, Gamma_N, fourier, no_gamma, U, c, t):
            
        u_ind, v_ind = self.V_ind_ub_field(x1_N, y1_N, x2_N, y2_N, Gamma_N, 0)

        # Finding induced velocity at each vortex by bound
        for n in range(len(u_ind)):
            
            u_ind_p, v_ind_p = self.V_ind_b_fast_4(fourier, x1_N[n], y1_N[n], c, t)

            u_ind[n] += u_ind_p #+ U_ref
            
            v_ind[n] += v_ind_p #+ pot.hdot(t) 


        return u_ind, v_ind
    
    def V_ind_ub_field(self, x1_N, y1_N, x2_N, y2_N, Gamma_N,vcore_flag):
        '''
        calculates induced velocity at (x1,y1) by vortices at (x2,y2)
        '''
        # Reshape the input arrays to enable broadcasting
        x1_N = np.reshape(x1_N, (-1, 1))  # Shape: (len(x1_N), 1)
        y1_N = np.reshape(y1_N, (-1, 1))  # Shape: (len(y1_N), 1)
        Gamma_N = np.reshape(Gamma_N, (-1, 1))  # Shape: (len(y1_N), 1)

        # Compute the difference arrays (x1_N - x2_N) and (y1_N - y2_N)
        dx = x1_N - x2_N  # Shape: (len(x1_N), len(x2_N))
        dy = y1_N - y2_N  # Shape: (len(y1_N), len(y2_N))

        # Calculate the squared distance and avoid division by zero by adding a small value (epsilon)
        r_squared = dx**2 + dy**2

        epsilon = 1e-10  # Small value to prevent division by zero
        r_squared = np.where(r_squared == 0, epsilon, r_squared)

        # Calculate induced velocities using broadcasting
        if vcore_flag == 1:
            u_ind =  dy*( 1 / (2 * np.pi * r_squared)) @ Gamma_N # Shape: (len(x1_N), len(x2_N))
            v_ind =  dx*(-1 / (2 * np.pi * r_squared)) @ Gamma_N     # Shape: (len(y1_N), len(y2_N))
        else:
            u_ind =  dy*( 1 / (2 * np.pi * (r_squared**2 + self.v_core**4)**0.5)) @ Gamma_N # Shape: (len(x1_N), len(x2_N))
            v_ind =  dx*(-1 / (2 * np.pi * (r_squared**2 + self.v_core**4)**0.5)) @ Gamma_N     # Shape: (len(y1_N), len(y2_N))

        return np.reshape(u_ind,-1), np.reshape(v_ind,-1)


    
    def fourier_gamma_calc_2(self, A_no, Gamma_N, eta_N, xi_N, N, c, t):

        fourier = np.zeros(A_no)

        theta = np.linspace(0,np.pi,200,endpoint=True)

        xi = 0.5*c*(1-np.cos(theta))

        eta = np.zeros(200)

        u_ind, v_ind = self.V_ind_ub_field(xi, eta, xi_N, eta_N, Gamma_N, 1)

        wx = - v_ind + self.kin.h_dot(t)

        fourier[0] = - 1 / np.pi / self.U_ref * inte.trapezoid(wx,theta)

        for i in range(1,A_no):

            fourier[i] = 2 / np.pi / self.U_ref * inte.trapezoid(wx*np.cos(i*theta),theta) 
            

        Gamma_b = np.pi * c * self.U_ref * (fourier[0] + fourier[1] * 0.5)
        # Gamma_b = self.U_ref * c * inte.trapezoid(fourier_inf,theta)

        return fourier, sum(Gamma_N), Gamma_b + sum(Gamma_N)

