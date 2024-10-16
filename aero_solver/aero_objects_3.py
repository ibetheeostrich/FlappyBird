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

        self.theta = linspace(0,pi,no_fourier*3,endpoint=True)
        self.x = 0.5*chords[0]*(1-cos(self.theta))
        self.y = zeros(no_fourier*3)

        self.uind = zeros(no_fourier*3)
        self.vind = zeros(no_fourier*3)
        self.dw   = zeros(no_fourier*3)

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

    def update_fourier_test(self, x_N, y_N,  Gamma_N, t):

        fourier = zeros(len(self.fourier))

        # v_core = 0.02*self.c(t)#1.3*self.t_step*self.x_dot(t)
        v_core = 1.3*self.t_step*self.x_dot(t)

        xi = 0.5 * self.c(t) * (1 - cos(self.theta))

        u_ind, v_ind = V_ind_ub_field(self.x, self.y, x_N, y_N, Gamma_N, v_core, 1)

        dphi_deta = v_ind*cos(self.alpha(t)) + u_ind*sin(self.alpha(t))

        wx = (- self.x_dot(t)*sin(self.alpha(t)) 
              - self.alpha_dot(t)*xi 
              + self.h_dot(t)*cos(self.alpha(t)) 
              - dphi_deta)

        fourier[0] = inte.trapezoid(- 1 / pi / self.x_dot(t) * wx,self.theta) 

        for i in range(1,len(self.fourier)):

            fourier[i] = inte.trapezoid(2 / pi / self.x_dot(t) * wx*cos(i*self.theta),self.theta)

        return fourier# - self.alpha(0.0)*self.x_dot(0.0)*pi*self.c(0.0)
    
    def update_fourier_1new(self, v_field, t):

        v_core = 1.3*self.t_step*self.x_dot(t)

        u_ind, v_ind = V_ind_ub_field(self.x, 
                                      self.y, 
                                      [v_field.tev_x[-1]], 
                                      [v_field.tev_y[-1]], 
                                      [v_field.tev[-1]],
                                      v_core, 1)

        dphi_deta = v_ind*cos(self.alpha(t)) + u_ind*sin(self.alpha(t))

        wx =  self.dw - dphi_deta

        self.fourier[0] = inte.trapezoid(- 1 / pi / self.x_dot(t) * wx ,self.theta) 

        for i in range(1,len(self.fourier)):

            self.fourier[i] = inte.trapezoid(2 / pi / self.x_dot(t) * wx*cos(i*self.theta),self.theta)

    
    def update_fourier_2new(self, v_field, t):

        v_core = 1.3*self.t_step*self.x_dot(t)

        u_ind, v_ind = V_ind_ub_field(self.x, 
                                      self.y, 
                                      [v_field.tev_x[-1], v_field.lev_x[-1]], 
                                      [v_field.tev_y[-1], v_field.lev_y[-1]], 
                                      [v_field.tev[-1], v_field.lev[-1]], 
                                      v_core, 1)

        dphi_deta = v_ind*cos(self.alpha(t)) + u_ind*sin(self.alpha(t))

        wx =  self.dw - dphi_deta

        self.fourier[0] = inte.trapezoid(- 1 / pi / self.x_dot(t) * wx ,self.theta) 

        for i in range(1,len(self.fourier)):

            self.fourier[i] = inte.trapezoid(2 / pi / self.x_dot(t) * wx * cos(i*self.theta),self.theta)

    
 
    def update_ind_vel(self, field, t):

        v_core = 1.3*self.t_step*self.x_dot(t)

        x1_N = self.x
        y1_N = self.y

        x2_N = concatenate((field.tev_x,field.lev_x,field.ext_x))
        y2_N = concatenate((field.tev_y,field.lev_y,field.ext_y))
        Gamma_N = concatenate((field.tev,field.lev,field.ext))

        self.uind, self.v_ind = V_ind_ub_field(x1_N, y1_N, x2_N, y2_N, Gamma_N, v_core,1)

    def update_downwash(self,t):

        xi = 0.5 * self.c(t) * (1 - cos(self.theta))

        dphi_deta = self.vind*cos(self.alpha(t)) + self.uind*sin(self.alpha(t))

        self.dw = (- self.x_dot(t)*sin(self.alpha(t)) 
                   - self.alpha_dot(t)*xi 
                   + self.h_dot(t)*cos(self.alpha(t)) 
                   - dphi_deta)

    def kelvinkutta_iter(self,inp ,v_field,t):

        # v_core = 0.02*self.c(t)#1.3*self.t_step*self.x_dot(t)
        v_core = 1.3*self.t_step*self.x_dot(t)

        u_ind, v_ind = V_ind_ub_field(self.x, 
                                      self.y, 
                                      [v_field.tev_x[-1]], 
                                      [v_field.tev_y[-1]], 
                                      [inp],
                                      v_core, 1)
        


        dphi_deta = v_ind*cos(self.alpha(t)) + u_ind*sin(self.alpha(t))

        wx = - dphi_deta

        a0 = inte.trapezoid(- 1 / pi / self.x_dot(t) * (wx + self.dw) ,self.theta)
        a1 = inte.trapezoid(  2 / pi / self.x_dot(t) * (wx + self.dw) * cos(self.theta),self.theta)
        gamma_b = pi * self.c(t) * self.x_dot(t) * (a0 + a1 * 0.5)
        gamma_b_old = pi * self.c(t) * self.x_dot(t) * (self.fourier_old[0] + self.fourier_old[1] * 0.5)

        gamma_field = sum(concatenate((v_field.tev[:-1],v_field.lev,v_field.ext)))

        return gamma_b - gamma_b_old + inp

    def kelvinkutta(self, v_field, t):

        iter = root(self.kelvinkutta_iter, 0.0, (v_field,t), tol=1e-16)

        v_field.tev[-1] = iter.x

        self.update_fourier_1new(v_field,t)


    def kelvin_lesp_iter(self,inp,v_field,lesp,t):

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

        a0 = inte.trapezoid(- 1 / pi / self.x_dot(t) * (wx + self.dw) ,self.theta)
        a1 = inte.trapezoid(  2 / pi / self.x_dot(t) * (wx + self.dw) * cos(self.theta),self.theta)
        gamma_b = pi * self.c(t) * self.x_dot(t) * (a0 + a1 * 0.5)
        gamma_b_old = pi * self.c(t) * self.x_dot(t) * (self.fourier_old[0] + self.fourier_old[1] * 0.5)


        gamma_field = sum(concatenate((v_field.tev[:-1],v_field.lev[:-1],v_field.ext)))

        return array([a0 - lesp, gamma_b + gamma_field + tev + lev])

    def kelvin_lesp_2(self,v_field,lesp,t):

        if lesp < 0:

            iter = root(self.kelvin_lesp_iter, array([1, -1]), (v_field,-lesp,t), tol=1e-16)
        else:
            iter = root(self.kelvin_lesp_iter, array([-1, 1]), (v_field,lesp,t), tol=1e-16)

        [v_field.tev[-1], v_field.lev[-1]] = iter.x
        

        self.update_fourier_2new(v_field, t)


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


        cnc  = 2.0*pi / self.x_dot(t) * (self.x_dot(t) * cos(self.alpha(t)) + 
                                           self.h_dot(t) * sin(self.alpha(t))) * (
                                               self.fourier[0] + 0.5 * self.fourier[1]
                                           )

        cnnc = 2.0*pi / self.x_dot(t)  * self.c(t) * (
            0.75  * (self.fourier[0] - self.fourier_old[0])/t_step + 
            0.25  * (self.fourier[1] - self.fourier_old[1])/t_step + 
            0.125 * (self.fourier[2] - self.fourier_old[2])/t_step
        )

        non1 = inte.trapezoid(2/self.x_dot(t)/self.x_dot(t)/self.c(t) * dphi_dxi*fourier_inf,self.theta)

        cn =  cnc + cnnc + non1

        # if t > 0.20:
        #     print(f"{cnc:.8f}", f"{cnnc:.8f}", f"{non1:.8f}", self.fourier[0], self.fourier[1])

        cs = 2*pi*self.fourier[0]**2

        cl = (0.5*1.225*self.x_dot(t)**2)*(cn*cos(self.alpha(t)) + cs*sin(self.alpha(t))) *self.c(t)
        cd = (0.5*1.225*self.x_dot(t)**2)*(-cn*sin(self.alpha(t)) + cs*cos(self.alpha(t)))*self.c(t)

        cl = (cn*cos(self.alpha(t)) + cs*sin(self.alpha(t))) 
        cd = (-cn*sin(self.alpha(t)) + cs*cos(self.alpha(t)))

        return cl, cd 
        # return (self.fourier[1] - self.fourier_old[1])/t_step, self.fourier[1]
        # return cnnc, non1
        # return (0.5*1.225*self.x_dot(t)**2)*cn*self.c(t), (0.5*1.225*self.x_dot(t)**2)*cs*self.c(t)

class vorticity_field:

    def __init__(self,c):

        self.tev = array([])
        self.lev = array([])
        self.ext = array([])

        self.tev_x = array([])
        self.lev_x = array([])
        self.ext_x = array([])

        self.tev_y = array([])
        self.lev_y = array([])
        self.ext_y = array([])

    def shed_tev(self, camber_line,t):

        if len(self.tev) == 0:
            self.tev_x = append(
                self.tev_x,
                camber_line.x[-1] + camber_line.x_dot(t)*0.5*camber_line.t_step
            )

            self.tev_y = append(
                self.tev_y,
                camber_line.y[-1] - camber_line.h_dot(t)*0.5*camber_line.t_step
            )

            self.tev = append(
                self.tev,
                0.0
            )

        else:
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
                    0.0
                )


            elif camber_line.fourier[0] < 0:
                self.lev = append(
                    self.lev,
                    0.0
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

    x1_N = reshape(x1_N, (len(x1_N),1))
    y1_N = reshape(y1_N, (len(x1_N),1))
    Gamma_N = reshape(Gamma_N, (len(Gamma_N),1))

    dx = x1_N - x2_N
    dy = y1_N - y2_N

    rsq = where(dx**2 + dy**2 == 0, 1, dx**2 + dy**2)

    inv = 0.5/sqrt(rsq**2 + v_core**4)/pi

    u_ind = dy * inv @ Gamma_N
    v_ind =-dx * inv @ Gamma_N

    return reshape(u_ind,-1), reshape(v_ind,-1)

def V_ind_b_fast_4(camber_line, x_n, y_n, v_core,t):

    v_core = 1.3*camber_line.t_step*camber_line.x_dot(t)

    fourier_inf = camber_line.fourier[0]*(1+cos(camber_line.theta))

    for i in range(1,len(camber_line.fourier)):

        fourier_inf += camber_line.fourier[i]*sin(i*camber_line.theta)*sin(camber_line.theta)

    fourier_inf *= camber_line.x_dot(t) * camber_line.c(t)
    # fourier_inf *= camber_line.V(t) * camber_line.c(t)


    inv = 1/sqrt(((x_n - camber_line.x)**2 + 
                  (y_n - camber_line.y)**2)**2 + 
                   v_core**4)
 

    x_s  = (x_n - camber_line.x) * inv 
    y_s  = (y_n - camber_line.y) * inv
    
    u_ind_p = 0.5 * inte.trapezoid(y_s*fourier_inf,camber_line.theta) / pi 

    v_ind_p =-0.5 * inte.trapezoid(x_s*fourier_inf,camber_line.theta) / pi

    return u_ind_p, v_ind_p



