# import numpy as np

from copy import deepcopy

from numpy import sin, cos, zeros, pi, linspace, concatenate, array, reshape, append, sqrt, where
from numpy.linalg import inv as inverse
from numpy.linalg import norm


import scipy.integrate as inte
from scipy.optimize import minimize
from scipy.optimize import root


err = 1e-12

class camber_line:

    def __init__(self,chords,no_fourier,x_dot,h_dot,alpha_dot,u,h,alpha,t_step):
    
        self.c = lambda t: chords[round(t/t_step)]
        self.le_pos = 0

        self.theta = linspace(0,pi,no_fourier*5,endpoint=True)
        self.x = 0.5*chords[0]*(1-cos(self.theta))
        self.y = zeros(no_fourier*5)

        self.fourier = zeros(no_fourier)
        self.fourier_old = zeros(no_fourier)

        self.x_dot = x_dot
        self.h_dot = h_dot
        self.alpha_dot = alpha_dot

        self.u = u
        self.h = h
        self.alpha = alpha

        self.V = lambda t: sqrt(x_dot(t)**2 + h_dot(t)**2)

        self.t_step = t_step

    def update_fourier(self, x_N, y_N,  Gamma_N, t):

        # v_core = 0.02*self.c(t)#1.3*self.t_step*self.x_dot(t)
        v_core = 1.3*self.t_step*self.x_dot(t)

        xi = 0.5 * self.c(t) * (1 - cos(self.theta))

        u_ind, v_ind = V_ind_ub_field(self.x, self.y, x_N, y_N, Gamma_N, v_core, 1)

        dphi_deta = v_ind*cos(self.alpha(t)) + u_ind*sin(self.alpha(t))

        wx = (- self.x_dot(t)*sin(self.alpha(t)) 
              - self.alpha_dot(t)*xi 
              + self.h_dot(t)*cos(self.alpha(t)) 
              - dphi_deta)

        self.fourier[0] = inte.trapezoid(- 1 / pi / self.x_dot(t) * wx,self.theta) 
        # self.fourier[0] = inte.trapezoid(- 1 / pi / self.V(t) * wx,self.theta) 


        for i in range(1,len(self.fourier)):

            self.fourier[i] = inte.trapezoid(2 / pi / self.x_dot(t) * wx*cos(i*self.theta),self.theta)
            # self.fourier[i] = inte.trapezoid(2 / pi / self.V(t) * wx*cos(i*self.theta),self.theta)

            

        Gamma_b = pi * self.c(t) * self.x_dot(t) * (self.fourier[0] + self.fourier[1] * 0.5)
        # Gamma_b = pi * self.c(t) * self.V(t) * (self.fourier[0] + self.fourier[1] * 0.5)

        # Gamma_b = self.U_ref * c * inte.trapezoid(fourier_inf,theta)

        return Gamma_b + sum(Gamma_N)# - self.alpha(0.0)*self.x_dot(0.0)*pi*self.c(0.0)
    
    def kelvinkutta_a0_a1(self, v_field, dh, t):

        # v_core = 0.02*self.c(t)#1.3*self.t_step*self.x_dot(t)
        v_core = 1.3*self.t_step*self.x_dot(t)

        u_ind, v_ind = V_ind_ub_field(self.x, 
                                      self.y, 
                                      v_field.tev_x[-1], 
                                      v_field.tev_y[-1], 
                                      v_field.tev[-1] + dh, 
                                      v_core, 1)

        dphi_deta = v_ind*cos(self.alpha(t)) + u_ind*sin(self.alpha(t))

        wx = - dphi_deta

        a0 = inte.trapezoid(- 1 / pi / self.x_dot(t) * wx,self.theta)
        a1 = inte.trapezoid(2 / pi / self.x_dot(t) * wx*cos(self.theta),self.theta)
        gamma_b = pi * self.c(t) * self.x_dot(t) * (a0 + a1 * 0.5)

        return gamma_b + v_field.tev[-1] + dh

    def kelvinkutta(self, v_field, dh, t):

        g0 = 100000

        while abs(g0) > err:

            g0 = self.update_fourier(concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                                     concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                                     concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                                     t)
            
            gp = self.kelvinkutta_a0_a1(v_field, dh, t)
            gm = self.kelvinkutta_a0_a1(v_field, -dh, t)
            
            v_field.tev[-1] = v_field.tev[-1] - dh * g0 / gp

            g0 = self.update_fourier(concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                                     concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                                     concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                                     t)


    def kelvinlesp_a0_a1(self, v_field, dg1, dg2, t):

        # v_core = 0.02*self.c(t)#1.3*self.t_step*self.x_dot(t)
        v_core = 1.3*self.t_step*self.x_dot(t)

        u_ind, v_ind = V_ind_ub_field(self.x, 
                                      self.y, 
                                      [v_field.tev_x[-1], v_field.lev_x[-1]], 
                                      [v_field.tev_y[-1], v_field.lev_y[-1]], 
                                      [v_field.tev[-1] + dg2, v_field.lev[-1] + dg1], 
                                      v_core, 1)

        dphi_deta = v_ind*cos(self.alpha(t)) + u_ind*sin(self.alpha(t))

        wx = - dphi_deta

        a0 = - 1 / pi / self.x_dot(t) * inte.trapezoid(wx,self.theta)
        a1 = 2 / pi / self.x_dot(t) * inte.trapezoid(wx*cos(self.theta),self.theta)
        gamma_b = pi * self.c(t) * self.x_dot(t) * (a0 + a1 * 0.5) 

        return gamma_b + v_field.tev[-1] + v_field.lev[-1] + dg1 + dg2, a0

    def kelvinlesp(self, v_field, dh, lesp, t):

        g0 = 100000

        if self.fourier[0] < 0:
            lesp_c = -lesp

        else:
            lesp_c = lesp

        iter = 1
        while abs(g0) > 1e-14 or abs(abs(self.fourier[0]) - lesp) > err:

            g0 = self.update_fourier(concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                                     concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                                     concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                                     t)

            g0_LEV_p, a0_LEV_p = self.kelvinlesp_a0_a1(v_field, dh, 0, t)          
            g0_LEV_m, a0_LEV_m = self.kelvinlesp_a0_a1(v_field,-dh, 0, t) 

            g0_TEV_p, a0_TEV_p = self.kelvinlesp_a0_a1(v_field, 0, dh, t)          
            g0_TEV_m, a0_TEV_m = self.kelvinlesp_a0_a1(v_field, 0,-dh, t)   

            target = array([self.fourier[0] - lesp_c, g0])

            jacob = array([[a0_LEV_p / dh,         
                               a0_TEV_p / dh],
                              [g0_LEV_p / dh, 
                               g0_TEV_p / dh]])
            
            try:

                # J_inv = inverse(jacob)

                # [v_field.lev[-1], v_field.tev[-1]] = array([v_field.lev[-1], v_field.tev[-1]]) - J_inv@target 

                # g0 = self.update_fourier(concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                #                          concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                #                          concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                #                          t)
                
                v_field.lev[-1] = v_field.lev[-1] + g0 * dh / g0_LEV_p - (self.fourier[0] + lesp_c) * dh / a0_LEV_p
                v_field.tev[-1] = v_field.tev[-1] + g0 * dh / g0_TEV_p - (self.fourier[0] + lesp_c) * dh / a0_TEV_p

                iter += 1
                if iter > 1000:
                    print("iter limit", g0, self.fourier[0])
                    break

            except:
                print('you are ugly and gay')
                return

    def kelvin_lesp_iter(self,inp,a0_int,g0_int,v_field,lesp,t):

        tev = inp[0]
        lev = inp[1]

        v_core = 1.3*self.t_step*self.x_dot(t)

        u_ind, v_ind = V_ind_ub_field(self.x, 
                                      self.y, 
                                      [v_field.tev_x[-1], v_field.lev_x[-1]], 
                                      [v_field.tev_y[-1], v_field.lev_y[-1]], 
                                      [tev, lev], 
                                      v_core, 1)
        
        dphi_deta = v_ind*cos(self.alpha(t)) + u_ind*sin(self.alpha(t))

        wx = - dphi_deta

        a0 = - 1 / pi / self.x_dot(t) * inte.trapezoid(wx,self.theta)
        a1 = 2 / pi / self.x_dot(t) * inte.trapezoid(wx*cos(self.theta),self.theta)
        gamma_b = pi * self.c(t) * self.x_dot(t) * (a0 + a1 * 0.5) 

        return [abs(gamma_b + tev + lev + g0_int), abs((a0 + a0_int) - lesp)]

    def kelvin_lesp_2(self,v_field,lesp,t):

        g0_int = deepcopy(self.update_fourier(concatenate((v_field.tev_x[:-1], v_field.lev_x[:-1], v_field.ext_x[:-1])),
                                     concatenate((v_field.tev_y[:-1], v_field.lev_y[:-1], v_field.ext_y[:-1])),
                                     concatenate((v_field.tev[:-1],   v_field.lev[:-1],   v_field.ext[:-1])),
                                     t))
            
        a0_int = deepcopy(self.fourier[0])


        iter = root(self.kelvin_lesp_iter, array([v_field.tev[-1], v_field.lev[-1]]), (a0_int,g0_int,v_field,lesp,t), tol=1e-12)

        [v_field.tev[-1], v_field.lev[-1]] = iter.x

        self.update_fourier(concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                            concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                            concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                            t)

    def update_pos(self,t):

        c = 0.5*self.c(t)*(1-cos(self.theta))

        self.x = cos(self.alpha(t)) * c - self.u(t)
        self.y =-sin(self.alpha(t)) * c + self.h(t)

    def calc_cl(self, x_N, y_N,  Gamma_N, t, t_step):

        # v_core = 0.02*self.c(t)#1.3*self.t_step*self.x_dot(t)
        v_core = 1.3*self.t_step*self.x_dot(t)


        u_ind, v_ind = V_ind_ub_field(self.x, self.y, x_N, y_N, Gamma_N, v_core, 1)

        dphi_dxi = -v_ind*sin(self.alpha(t)) + u_ind*cos(self.alpha(t))

        fourier_inf = self.fourier[0]*(1+cos(self.theta))

        for i in range(1,len(self.fourier)):

            fourier_inf += self.fourier[i]*sin(i*self.theta)*sin(self.theta)

        fourier_inf *= self.x_dot(t) * self.c(t)
        # fourier_inf *= self.V(t) * self.c(t)


        cnc = 2.0*pi / self.x_dot(t) * (self.x_dot(t) * cos(self.alpha(t)) + 
                                           self.h_dot(t) * sin(self.alpha(t))) * (
                                               self.fourier[0] + 0.5 * self.fourier[1]
                                           )
        
        # cnnc = 2.0*pi * self.c(t) * (
        #     0.75  * (self.x_dot(t) * self.fourier[0] - self.x_dot(t-t_step) * self.fourier_old[0])/t_step + 
        #     0.25  * (self.x_dot(t) * self.fourier[1] - self.x_dot(t-t_step) * self.fourier_old[1])/t_step + 
        #     0.125 * (self.x_dot(t) * self.fourier[2] - self.x_dot(t-t_step) * self.fourier_old[2])/t_step
        # )

        cnnc = 2.0*pi * self.c(t) * (
            0.75  * (self.fourier[0] - self.fourier_old[0])/t_step + 
            0.25  * (self.fourier[1] - self.fourier_old[1])/t_step + 
            0.125 * (self.fourier[2] - self.fourier_old[2])/t_step
        )

        non1 = inte.trapezoid(2/self.x_dot(t)/self.x_dot(t)/self.c(t) * dphi_dxi*fourier_inf,self.theta)

        cn = cnc + cnnc + non1

        # if t > 0.20:
        #     print(f"{cnc:.8f}", f"{cnnc:.8f}", f"{non1:.8f}", self.fourier[0], self.fourier[1])

        cs = 2*pi*self.fourier[0]**2

        cl = (0.5*1.225*self.x_dot(t)**2)*(cn*cos(self.alpha(t)) + cs*sin(self.alpha(t))) *self.c(t)
        cd = (0.5*1.225*self.x_dot(t)**2)*(-cn*sin(self.alpha(t)) + cs*cos(self.alpha(t)))*self.c(t)

        return cl, cd 
        # return (self.fourier[1] - self.fourier_old[1])/t_step, self.fourier[1]
        # return cnnc, non1
        # return (0.5*1.225*self.x_dot(t)**2)*cn*self.c(t), (0.5*1.225*self.x_dot(t)**2)*cs*self.c(t)

class vorticity_field:

    def __init__(self,c):

        self.tev = array([0.0])
        self.lev = array([])
        self.ext = array([])

        self.tev_x = array([c])
        self.lev_x = array([])
        self.ext_x = array([])

        self.tev_y = array([0.0])
        self.lev_y = array([])
        self.ext_y = array([])

    def shed_tev(self, camber_line):

        self.tev_x = append(
            self.tev_x,
            camber_line.x[-1] + 
            (self.tev_x[-1] - camber_line.x[-1])/3
        )
        
        self.tev_y = append(
            self.tev_y,
            camber_line.y[-1] + 
            (self.tev_y[-1] - camber_line.y[-1])/3
        )
        
        self.tev = append(
            self.tev,
            0.0
        )

    def shed_lev(self, camber_line,t):

        if len(self.lev) == 0:
        # if True:
            self.lev_x = append(
                self.lev_x,
                camber_line.x[0] +
                (camber_line.x[0] - camber_line.x[1])*0.05
            )

            self.lev_y = append(
                self.lev_y,
                camber_line.y[0]+
                (camber_line.y[0] - camber_line.y[1])*0.05
            )

            # self.lev_x = append(
            #     self.lev_x,
            #     camber_line.x[0] + camber_line.x_dot(t)*0.5*camber_line.t_step
            # )

            # self.lev_y = append(
            #     self.lev_y,
            #     camber_line.y[0] - camber_line.h_dot(t)*0.5*camber_line.t_step
            # )

            if camber_line.fourier[0] > 0:
                self.lev = append(
                    self.lev,
                    10.0
                )


            elif camber_line.fourier[0] < 0:
                self.lev = append(
                    self.lev,
                    -10.0
                )

        else:

            self.lev_x = append(
                self.lev_x,
                camber_line.x[0] +
                (self.lev_x[-1] - camber_line.x[0])/3
            )

            self.lev_y = append(
                self.lev_y,
                camber_line.y[0] +
                (self.lev_y[-1] - camber_line.y[0])/3
            
            
            )
            if camber_line.fourier[0] > 0:
                self.lev = append(
                    self.lev,
                    10
                )

            elif camber_line.fourier[0] < 0:
                self.lev = append(
                    self.lev,
                    -10
                )

    def advect(self, camber_line,t_step,t):

        # v_core = 0.02*camber_line.c(t)
        v_core = 1.3*t_step*camber_line.x_dot(t)

        x_tot = concatenate((self.tev_x, self.lev_x, self.ext_x))
        y_tot = concatenate((self.tev_y, self.lev_y, self.ext_y))
        g_tot = concatenate((self.tev,   self.lev,   self.ext))

        u_ind, v_ind = V_ind_ub_field(
            x_tot,
            y_tot,
            x_tot,
            y_tot,
            g_tot,
            v_core,
            1)

        for i in range(len(u_ind)):

            u_ind_p, v_ind_p = V_ind_b_fast_4(camber_line, x_tot[i], y_tot[i], v_core,t)
            u_ind[i] += u_ind_p
            v_ind[i] += v_ind_p

        x_tot += u_ind*t_step
        y_tot += v_ind*t_step

        self.tev_x = x_tot[0:len(self.tev_x)]
        self.tev_y = y_tot[0:len(self.tev_x)]
        self.tev   = g_tot[0:len(self.tev_x)]

        self.lev_x = x_tot[len(self.tev_x) : len(self.tev_x) + len(self.lev_x)]
        self.lev_y = y_tot[len(self.tev_x) : len(self.tev_x) + len(self.lev_x)]
        self.lev   = g_tot[len(self.tev_x) : len(self.tev_x) + len(self.lev_x)]

        self.ext_x = x_tot[len(self.tev_x) + len(self.lev_x) : len(self.tev_x) + len(self.lev_x) + len(self.ext_x)]
        self.ext_y = y_tot[len(self.tev_x) + len(self.lev_x) : len(self.tev_x) + len(self.lev_x) + len(self.ext_x)]
        self.ext   = g_tot[len(self.tev_x) + len(self.lev_x) : len(self.tev_x) + len(self.lev_x) + len(self.ext_x)]

def V_ind_ub_field(x1_N, y1_N, x2_N, y2_N, Gamma_N, v_core,v_core_flag):
    '''
    calculates induced velocity at (x1,y1) by vortices at (x2,y2)
    '''
    # Reshape the input arrays to enable broadcasting
    x1_N = reshape(x1_N, (-1, 1))  # Shape: (len(x1_N), 1)
    y1_N = reshape(y1_N, (-1, 1))  # Shape: (len(y1_N), 1)
    Gamma_N = reshape(Gamma_N, (-1, 1))  # Shape: (len(y1_N), 1)

    # Compute the difference arrays (x1_N - x2_N) and (y1_N - y2_N)
    dx = x1_N - x2_N  # Shape: (len(x1_N), len(x2_N))
    dy = y1_N - y2_N  # Shape: (len(y1_N), len(y2_N))

    # Calculate the squared distance and avoid division by zero by adding a small value (epsilon)
    r_squared = dx**2 + dy**2

    epsilon = 1e-10  # Small value to prevent division by zero
    r_squared = where(r_squared == 0, epsilon, r_squared)

    # Calculate induced velocities using broadcasting
    if v_core_flag == 0:
        u_ind =  dy*( 1 / (2 * pi * r_squared)) @ Gamma_N # Shape: (len(x1_N), len(x2_N))
        v_ind =  dx*(-1 / (2 * pi * r_squared)) @ Gamma_N     # Shape: (len(y1_N), len(y2_N))
    else:
        u_ind =  dy*( 1 / (2 * pi * sqrt(r_squared**2 + v_core**4))) @ Gamma_N # Shape: (len(x1_N), len(x2_N))
        v_ind =  dx*(-1 / (2 * pi * sqrt(r_squared**2 + v_core**4))) @ Gamma_N     # Shape: (len(y1_N), len(y2_N))

    return reshape(u_ind,-1), reshape(v_ind,-1)

def V_ind_b_fast_4(camber_line, x_n, y_n, v_core,t):

    v_core = 1.3*camber_line.t_step*camber_line.x_dot(t)

    fourier_inf = camber_line.fourier[0]*(1+cos(camber_line.theta))

    for i in range(1,len(camber_line.fourier)):

        fourier_inf += camber_line.fourier[i]*sin(i*camber_line.theta)*sin(camber_line.theta)

    # fourier_inf *= camber_line.x_dot(t) * camber_line.c(t)
    fourier_inf *= camber_line.V(t) * camber_line.c(t)


    inv = 1/sqrt(((x_n - camber_line.x)**2 + 
                  (y_n - camber_line.y)**2)**2 + 
                   v_core**4)
 

    x_s  = (x_n - camber_line.x) * inv
    y_s  = (y_n - camber_line.y) * inv
    
    u_ind = 0.5 * inte.trapezoid(y_s*fourier_inf,camber_line.theta) / pi

    v_ind =-0.5 * inte.trapezoid(x_s*fourier_inf,camber_line.theta) / pi

    return u_ind, v_ind



