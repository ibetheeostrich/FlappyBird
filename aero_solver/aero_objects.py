import math as mt
import numpy as np

import scipy.integrate as inte
import matplotlib.pyplot as plt


class camber_line:

    def __init__(self,chords,no_fourier,x_dot,h_dot,alpha_dot,u,h,alpha,t_step):
    
        self.c = lambda t: chords[round(t/t_step)]
        self.le_pos = 0
        self.alpha = 0

        self.theta = np.linspace(0,np.pi,no_fourier*2,endpoint=True)
        self.x = 0.5*chords[0]*(1-np.cos(self.theta))
        self.y = np.zeros(no_fourier*2)

        self.fourier = np.zeros(no_fourier)

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

        return sum(Gamma_N), Gamma_b + sum(Gamma_N)
    
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

            g0 = self.kelvinkutta_a0_a1(v_field, 0, t)
            gp = self.kelvinkutta_a0_a1(v_field, dh, t)
            gm = self.kelvinkutta_a0_a1(v_field, -dh, t)
            
            v_field.tev[-1] = v_field.tev[-1] - 2*dh * g0 / (gp-gm)

        

    def update_pos(self,t):

        c = 0.5*self.c(t)*(1-np.cos(self.theta))

        self.x = np.cos(self.alpha(t)) * c - self.u(t)
        self.y =-np.sin(self.alpha(t)) * c + self.h(t)

    def calc_cl(self):

        return np.pi*(2*self.fourier[0] + self.fourier[1])

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

    def advect(self, camber_line,t_step):

        v_core = 0.02 * camber_line.c(t)

        x_tot = np.concatenate((field.tev_x, field.lev_x, field.ext_x))
        y_tot = np.concatenate((field.tev_y, field.lev_y, field.ext_y))
        g_tot = np.concatenate((field.tev, field.lev, field.ext))

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

        field.tev_x = x_tot[0:len(field.tev_x)]
        field.tev_y = y_tot[0:len(field.tev_x)]
        field.tev   = g_tot[0:len(field.tev_x)]

        field.lev_x = x_tot[len(field.tev_x) : len(field.tev_x) + len(field.lev_x)]
        field.lev_y = y_tot[len(field.tev_x) : len(field.tev_x) + len(field.lev_x)]
        field.lev   = g_tot[len(field.tev_x) : len(field.tev_x) + len(field.lev_x)]

        field.ext_x = x_tot[len(field.tev_x) + len(field.lev_x) : len(field.tev_x) + len(field.lev_x) + len(field.ext_x)]
        field.ext_y = y_tot[len(field.tev_x) + len(field.lev_x) : len(field.tev_x) + len(field.lev_x) + len(field.ext_x)]
        field.ext   = g_tot[len(field.tev_x) + len(field.lev_x) : len(field.tev_x) + len(field.lev_x) + len(field.ext_x)]

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

chords = 1 + np.zeros(400)
t_step = 0.01
td = np.linspace(0,400*t_step,400,endpoint=False)


x_dot = lambda t: 5
h_dot = lambda t: 0
alpha_dot = lambda t: 2*np.pi*0.5*np.sin(2*np.pi*t)

u = lambda t: 5*t
h = lambda t: 0
alpha = lambda t: 0.5 - 0.5*np.cos(2*np.pi*t)

bem = camber_line(chords, 50, x_dot,h_dot,alpha_dot,u,h,alpha,t_step)

field = vorticity_field(chords[0])

for t in td:

    if t > 0:
        
        bem.update_pos(t)

        field.shed_tev(bem)

        bem.kelvinkutta(field,0.001,t)

        bem.update_fourier(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                           np.concatenate((field.tev_y, field.lev_y, field.ext_y)),
                           np.concatenate((field.tev, field.lev, field.ext)),
                           t)
        
        # field.advect(bem,t_step)

        print(bem.calc_cl())



fig, ax = plt.subplots()
fig.dpi = 300
ax.plot(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
        np.concatenate((field.tev_y, field.lev_y, field.ext_y))
        ,'ro')
ax.plot(bem.x,
        bem.y,
        'k')
ax.axis("equal")
plt.savefig('test1' + '.png')
plt.clf()