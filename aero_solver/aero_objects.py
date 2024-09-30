import math as mt
import numpy as np

import scipy.integrate as inte
import matplotlib.pyplot as plt

from copy import deepcopy


class camber_line:

    def __init__(self,chords,no_fourier,x_dot,h_dot,alpha_dot,u,h,alpha,t_step):
    
        self.c = lambda t: chords[round(t/t_step)]
        self.le_pos = 0
        self.alpha = 0

        self.theta = np.linspace(0,np.pi,no_fourier*2,endpoint=True)
        self.x = 0.5*chords[0]*(1-np.cos(self.theta))
        self.y = np.zeros(no_fourier*2)

        self.fourier = np.zeros(no_fourier)
        self.fourier_old = np.zeros(no_fourier)

        self.x_dot = x_dot
        self.h_dot = h_dot
        self.alpha_dot = alpha_dot

        self.u = u
        self.h = h
        self.alpha = alpha

    def update_fourier(self, x_N, y_N,  Gamma_N, t):

        v_core = 0.02 * self.c(t)
        xi = 0.5 * self.c(t) * (1 - np.cos(self.theta))

        u_ind, v_ind = V_ind_ub_field(self.x, self.y, x_N, y_N, Gamma_N, v_core, 1)

        dphi_deta = v_ind*np.cos(self.alpha(t)) + u_ind*np.sin(self.alpha(t))

        wx = (- self.x_dot(t)*np.sin(self.alpha(t)) 
              - self.alpha_dot(t)*xi 
              + self.h_dot(t)*np.cos(self.alpha(t)) 
              - dphi_deta)

        self.fourier[0] = - 1 / np.pi / self.x_dot(t) * inte.trapezoid(wx,self.theta)

        for i in range(1,len(self.fourier)):

            self.fourier[i] = 2 / np.pi / self.x_dot(t) * inte.trapezoid(wx*np.cos(i*self.theta),self.theta) 
            

        Gamma_b = np.pi * self.c(t) * self.x_dot(t) * (self.fourier[0] + self.fourier[1] * 0.5)
        # Gamma_b = self.U_ref * c * inte.trapezoid(fourier_inf,theta)

        return Gamma_b + sum(Gamma_N)
    
    def kelvinkutta_a0_a1(self, v_field, dh, t):

        v_core = 0.02 * self.c(t)

        u_ind, v_ind = V_ind_ub_field(self.x, 
                                      self.y, 
                                      v_field.tev_x[-1], 
                                      v_field.tev_y[-1], 
                                      v_field.tev[-1] + dh, 
                                      v_core, 1)

        dphi_deta = v_ind*np.cos(self.alpha(t)) + u_ind*np.sin(self.alpha(t))

        wx = - dphi_deta

        a0 = - 1 / np.pi / self.x_dot(t) * inte.trapezoid(wx,self.theta)
        a1 = 2 / np.pi / self.x_dot(t) * inte.trapezoid(wx*np.cos(self.theta),self.theta)
        gamma_b = np.pi * self.c(t) * self.x_dot(t) * (a0 + a1 * 0.5) 

        return gamma_b + v_field.tev[-1] + dh

    def kelvinkutta(self, v_field, dh, t):

        g0 = 100000

        while abs(g0) > 0.00001:

            g0 = self.update_fourier(np.concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                                     np.concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                                     np.concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                                     t)
            
            gp = self.kelvinkutta_a0_a1(v_field, dh, t)
            gm = self.kelvinkutta_a0_a1(v_field, -dh, t)
            
            v_field.tev[-1] = v_field.tev[-1] - 2*dh * g0 / (gp-gm)

            g0 = self.update_fourier(np.concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                                     np.concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                                     np.concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                                     t)


    def kelvinlesp_a0_a1(self, v_field, dg1, dg2, t):

        v_core = 0.02 * self.c(t)

        u_ind, v_ind = V_ind_ub_field(self.x, 
                                      self.y, 
                                      [v_field.tev_x[-1], v_field.lev_x[-1]], 
                                      [v_field.tev_y[-1], v_field.lev_y[-1]], 
                                      [v_field.tev[-1] + dg2, v_field.lev[-1] + dg1], 
                                      v_core, 1)

        dphi_deta = v_ind*np.cos(self.alpha(t)) + u_ind*np.sin(self.alpha(t))

        wx = - dphi_deta

        a0 = - 1 / np.pi / self.x_dot(t) * inte.trapezoid(wx,self.theta)
        a1 = 2 / np.pi / self.x_dot(t) * inte.trapezoid(wx*np.cos(self.theta),self.theta)
        gamma_b = np.pi * self.c(t) * self.x_dot(t) * (a0 + a1 * 0.5) 

        return gamma_b + v_field.tev[-1] + v_field.lev[-1] + dg1 + dg2, a0

    def kelvinlesp(self, v_field, dh, lesp, t):

        g0 = 100000

        if self.fourier[0] < 0:
            lesp_c = -lesp

        else:
            lesp_c = lesp


        while abs(g0) > 0.0000001 and abs(abs(self.fourier[0]) - lesp) > 0.00000001:

            g0 = self.update_fourier(np.concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                                     np.concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                                     np.concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                                     t)

            g0_LEV_p, a0_LEV_p = self.kelvinlesp_a0_a1(v_field, dh, 0, t)          
            g0_LEV_m, a0_LEV_m = self.kelvinlesp_a0_a1(v_field,-dh, 0, t) 

            g0_TEV_p, a0_TEV_p = self.kelvinlesp_a0_a1(v_field, 0, dh, t)          
            g0_TEV_m, a0_TEV_m = self.kelvinlesp_a0_a1(v_field, 0,-dh, t)   

            target = np.array([self.fourier[0] - lesp_c, g0])

            jacob = np.array([[(a0_LEV_p - a0_LEV_m) / (2*dh),         
                               (a0_TEV_p - a0_TEV_m) / (2*dh)],
                              [(g0_LEV_p - g0_LEV_m) / (2*dh), 
                               (g0_TEV_p - g0_TEV_m) / (2*dh)]])
            
            try:

                J_inv = np.linalg.inv(jacob)

                [v_field.lev[-1], v_field.tev[-1]] = np.array([v_field.lev[-1], v_field.tev[-1]]) - J_inv@target 

                g0 = self.update_fourier(np.concatenate((v_field.tev_x, v_field.lev_x, v_field.ext_x)),
                                         np.concatenate((v_field.tev_y, v_field.lev_y, v_field.ext_y)),
                                         np.concatenate((v_field.tev,   v_field.lev,   v_field.ext)),
                                         t)

            except:
                print('you are ugly and gay')
                return

    def update_pos(self,t):

        c = 0.5*self.c(t)*(1-np.cos(self.theta))

        self.x = np.cos(self.alpha(t)) * c - self.u(t)
        self.y =-np.sin(self.alpha(t)) * c + self.h(t)

    def calc_cl(self, x_N, y_N,  Gamma_N, t, t_step):

        v_core = 0.02 * self.c(t)
        xi = 0.5 * self.c(t) * (1 - np.cos(self.theta))

        u_ind, v_ind = V_ind_ub_field(self.x, self.y, x_N, y_N, Gamma_N, v_core, 1)

        dphi_dxi = -v_ind*np.sin(self.alpha(t)) + u_ind*np.cos(self.alpha(t))

        fourier_inf = self.fourier[0]*(1+np.cos(self.theta))

        for i in range(1,len(self.fourier)):

            fourier_inf += self.fourier[i]*np.sin(i*self.theta)*np.sin(self.theta)

        fourier_inf *= self.x_dot(t) * self.c(t)

        cnc = 2.0*np.pi / self.x_dot(t) * (self.x_dot(t) * np.cos(self.alpha(t)) + 
                                           self.h_dot(t)*np.sin(self.alpha(t))) * (
                                               self.fourier[0] + 0.5 * self.fourier[1]
                                           )
        
        cnnc = 2.0*np.pi / self.x_dot(t) / self.c(t) * (
            0.75  * (self.fourier[0] - self.fourier_old[0]) / t_step + 
            0.25  * (self.fourier[1] - self.fourier_old[1]) / t_step + 
            0.125 * (self.fourier[2] - self.fourier_old[2]) / t_step
        )

        non1 = 2/self.x_dot(t)/self.x_dot(t)/self.c(t) * inte.trapezoid(dphi_dxi*dphi_dxi,self.theta)

        cn = cnc + cnnc + non1

        cs = 2*np.pi*self.fourier[0]**2

        # print(cnc, cnnc, non1, self.alpha(t))

        return (self.c(t)*0.5*1.225*self.x_dot(t)**2)*(cn*np.cos(self.alpha(t)) + cs*np.sin(self.alpha(t)))

class vorticity_field:

    def __init__(self,c):

        self.tev = np.array([0.0])
        self.lev = np.array([])
        self.ext = np.array([])

        self.tev_x = np.array([c])
        self.lev_x = np.array([])
        self.ext_x = np.array([])

        self.tev_y = np.array([0.0])
        self.lev_y = np.array([])
        self.ext_y = np.array([])

    def shed_tev(self, camber_line):

        self.tev_x = np.append(
            self.tev_x,
            camber_line.x[-1] + 
            (self.tev_x[-1] - camber_line.x[-1])/3
        )
        
        self.tev_y = np.append(
            self.tev_y,
            camber_line.y[-1] + 
            (self.tev_y[-1] - camber_line.y[-1])/3
        )
        
        self.tev = np.append(
            self.tev,
            0.0
        )

    def shed_lev(self, camber_line):

        self.lev_x = np.append(
            self.lev_x,
            camber_line.x[0]
        )
    
        
        if camber_line.fourier[0] > 0:
            self.lev = np.append(
                self.lev,
                10
            )

            self.lev_y = np.append(
                self.lev_y,
                camber_line.y[0]
            )

        elif camber_line.fourier[0] < 0:
            self.lev = np.append(
                self.lev,
                -10
            )

            self.lev_y = np.append(
                self.lev_y,
                camber_line.y[0]
            )

    def advect(self, camber_line,t_step,t):

        v_core = 0.02 * camber_line.c(t)

        x_tot = np.concatenate((self.tev_x, self.lev_x, self.ext_x))
        y_tot = np.concatenate((self.tev_y, self.lev_y, self.ext_y))
        g_tot = np.concatenate((self.tev,   self.lev,   self.ext))

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
    if v_core_flag == 0:
        u_ind =  dy*( 1 / (2 * np.pi * r_squared)) @ Gamma_N # Shape: (len(x1_N), len(x2_N))
        v_ind =  dx*(-1 / (2 * np.pi * r_squared)) @ Gamma_N     # Shape: (len(y1_N), len(y2_N))
    else:
        u_ind =  dy*( 1 / (2 * np.pi * (r_squared**2 + v_core**4)**0.5)) @ Gamma_N # Shape: (len(x1_N), len(x2_N))
        v_ind =  dx*(-1 / (2 * np.pi * (r_squared**2 + v_core**4)**0.5)) @ Gamma_N     # Shape: (len(y1_N), len(y2_N))

    return np.reshape(u_ind,-1), np.reshape(v_ind,-1)

def V_ind_b_fast_4(camber_line, x_n, y_n, v_core,t):


    fourier_inf = camber_line.fourier[0]*(1+np.cos(camber_line.theta))

    for i in range(1,len(camber_line.fourier)):

        fourier_inf += camber_line.fourier[i]*np.sin(i*camber_line.theta)*np.sin(camber_line.theta)

    x_s  = (x_n - camber_line.x) / np.sqrt(((x_n - camber_line.x)**2 + 
                                            (y_n - camber_line.y)**2)**2 
                                            + v_core**4)

    y_s  = (y_n - camber_line.y) / np.sqrt(((x_n - camber_line.x)**2 + 
                                            (y_n - camber_line.y)**2)**2 + 
                                            v_core**4)
 
    
    u_ind = 0.5 * camber_line.c(t) * camber_line.x_dot(t) * inte.trapezoid(y_s*fourier_inf,camber_line.theta) / np.pi

    v_ind =-0.5 * camber_line.c(t) * camber_line.x_dot(t) * inte.trapezoid(x_s*fourier_inf,camber_line.theta) / np.pi

    return u_ind, v_ind



