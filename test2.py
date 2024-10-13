import math
import aero_solver.pot_func_simp as aero
import numpy as np
from copy import deepcopy
import scipy.integrate as inte
import aero_solver.aero_objects as ao

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import io


from aero_solver.bem import bem as bem

t_step = 0.05
no_steps = 800
chords = 1.0+np.zeros(no_steps)#*np.linspace(0,0.25,no_steps)
freq = 3
amp = np.deg2rad(42.5)


tag = 9
cl = np.zeros(no_steps)
cd = np.zeros(no_steps)

td = np.linspace(0,no_steps*t_step,no_steps,endpoint=False)

lesp_crit = 0.195

x_dot = lambda t: 1
h_dot = lambda t: 0.0#2*np.pi*np.sin(2*np.pi*t)
alpha_dot = lambda t: 0.0

u = lambda t: 1.0*t
h = lambda t: 0.0#1-np.cos(2*np.pi*t)
alpha = lambda t: np.deg2rad(10)

be  = ao.camber_line(chords, 35, x_dot,h_dot,alpha_dot,u,h,alpha,t_step)

field   = ao.vorticity_field(chords[0])



lev_flag = 0

for t in td:

    if t > -1:
        lev_flag_prev = lev_flag

        if t > t_step:
            be.fourier_old = deepcopy(be.fourier) 

        be.update_pos(t)

        field.shed_tev(be)

        be.kelvinkutta(field,0.001,t)

        if abs(be.fourier[0]) > lesp_crit:          

            field.shed_lev(be,t)

            # be.kelvinlesp(field, 0.001, lesp_crit, t)
            be.kelvin_lesp_2(field,lesp_crit,t)

            gb = np.pi * be.c(t) * be.x_dot(t) * (be.fourier[0] + be.fourier[1] * 0.5) + np.sum(np.concatenate((field.tev, field.lev, field.ext)))           
            # print(f'{t:.2f} {abs(be.fourier[0]) - lesp_crit:.3f} {gb}') 

            lev_flag = 1

        if lev_flag_prev == 1 and lev_flag == 0:

            field.ext   = np.append(field.ext   ,deepcopy(field.lev)    )
            field.ext_x = np.append(field.ext_x ,deepcopy(field.lev_x)  )
            field.ext_y = np.append(field.ext_y ,deepcopy(field.lev_y)  )

            field.lev   = np.array([])
            field.lev_x = np.array([])
            field.lev_y = np.array([])
#####################################################################################    

        if round(t/t_step) % 5== 0:
        
            fig, ax = plt.subplots()
            fig.dpi = 300
            ax.scatter(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                    np.concatenate((field.tev_y, field.lev_y, field.ext_y))
                    , c='b', s=1.7)
            ax.plot(be.x,
                    be.y,
                    'k')
            ax.axis("equal")
            # ax.set_xlim(be.x[0] - 0.1,be.x[-1] + 0.1)
            # ax.set_ylim(be.y[0] - 0.1,be.y[-1] + 0.1)
            # plt.savefig('./results/' + str(tag) + '/' + str(round(t/t_step)) + '.png')
            plt.savefig(str(round(t/t_step)) + '.png')
            plt.clf()   

            gb = np.pi * be.c(t) * be.x_dot(t) * (be.fourier[0] + be.fourier[1] * 0.5) + np.sum(np.concatenate((field.tev, field.lev, field.ext)))           
            # print(f'{t:.2f} {abs(be.fourier[0]) - lesp_crit:.3f} {gb}') 
        # print(x_dot(t))

#####################################################################################    

        field.advect(be,t_step,t)

        cl[round(t/t_step)], cd[round(t/t_step)] = be.calc_cl(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                                         np.concatenate((field.tev_y, field.lev_y, field.ext_y)),
                                         np.concatenate((field.tev, field.lev, field.ext)),
                                         t, t_step)

        print(cl[round(t/t_step)])


fig, ax = plt.subplots()
fig.dpi = 300
ax.plot(td[2:],
        cl[2:],
        'k')
plt.savefig('flat1' + '.png')
plt.clf()   

np.savetxt(str(10) + 'deg.csv', np.transpose(np.vstack((cl,td))),  
          delimiter = ",")