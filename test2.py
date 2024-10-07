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

t_step = 0.0015
no_steps = 800
chords = 1.0+np.zeros(no_steps)#*np.linspace(0,0.25,no_steps)
freq = 3
amp = np.deg2rad(42.5)


tag = 9
cl = np.zeros(no_steps)
cd = np.zeros(no_steps)

td = np.linspace(0,no_steps*t_step,no_steps,endpoint=False)

lesp_crit = 0.19

x_dot = lambda t: 5.0
h_dot = lambda t: 2*np.pi*np.sin(2*np.pi*t)
alpha_dot = lambda t: 0.0

u = lambda t: 5.0*t
h = lambda t: 1-np.cos(2*np.pi*t)
alpha = lambda t: 0.0

be  = ao.camber_line(chords, 35, x_dot,h_dot,alpha_dot,u,h,alpha,t_step)

field   = ao.vorticity_field(chords[0])

lev_flag = 0

for t in td:

    if t > 0.0:

        lev_flag_prev = lev_flag

        be.fourier_old = deepcopy(be.fourier) 

        be.update_pos(t)

        field.shed_tev(be)

        be.kelvinkutta(field,0.000001,t)

        lev_flag = 0

        if abs(be.fourier[0]) > lesp_crit:          

            field.shed_lev(be)

            be.kelvinlesp(field, 0.000001, lesp_crit, t)

            lev_flag = 1

        if lev_flag_prev == 1 and lev_flag == 0:

            field.ext   = np.append(field.ext   ,deepcopy(field.lev)    )
            field.ext_x = np.append(field.ext_x ,deepcopy(field.lev_x)  )
            field.ext_y = np.append(field.ext_y ,deepcopy(field.lev_y)  )

            field.lev   = np.array([])
            field.lev_x = np.array([])
            field.lev_y = np.array([])

#####################################################################################    

        # if round(t/t_step) % 5 == 0 and tag == 9:
        #     fig, ax = plt.subplots()
        #     fig.dpi = 300
        #     ax.scatter(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
        #             np.concatenate((field.tev_y, field.lev_y, field.ext_y))
        #             , c='b', s=1.7)
        #     ax.plot(be.x,
        #             be.y,
        #             'k')
        #     ax.axis("equal")
        #     # ax.set_xlim(be.x[0] - 0.1,be.x[-1] + 0.1)
        #     # ax.set_ylim(be.y[0] - 0.1,be.y[-1] + 0.1)
        #     plt.savefig(str(round(t/t_step)) + '.png')
        #     plt.clf()   

                        
        t1 = np.pi * be.c(t) * be.x_dot(t) * (be.fourier[0] + be.fourier[1] * 0.5) + np.sum(np.concatenate((field.tev, field.lev, field.ext)))
        t2 = be.fourier[0]

        print(t1,t2)

#####################################################################################    

        field.advect(be,t_step,t)

        cl[round(t/t_step)], cd[round(t/t_step)] = be.calc_cl(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                                         np.concatenate((field.tev_y, field.lev_y, field.ext_y)),
                                         np.concatenate((field.tev, field.lev, field.ext)),
                                         t, t_step)

# Creating Figures
fig, ax = plt.subplots()
fig.dpi = 300
ax.plot(td, cl)
ax.plot(td, cd)
ax.set_xlabel('Time  (s)')
ax.set_ylabel('Lift Force (n)')
plt.savefig('long_test' + '.png')
plt.close(fig)
