import math
import aero_solver.pot_func_simp as aero
import numpy as np
from copy import deepcopy
import scipy.integrate as inte
import aero_solver.aero_objects_2 as ao

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import io


from aero_solver.bem import bem_2 as bem

csfont = {'fontname':'Comic Sans MS'}
hfont = {'fontname':'Helvetica'}

t_step = 0.005
no_steps = 800
chords = 1.0+np.zeros(no_steps)#*np.linspace(0,0.25,no_steps)
freq = 3
amp = np.deg2rad(42.5)


tag = 9
cl = np.zeros(no_steps)
cd = np.zeros(no_steps)

td = np.linspace(0,no_steps*t_step,no_steps,endpoint=False)

lesp_crit = 0.2

x_dot = lambda t: 5
h_dot = lambda t: 0.0#2*np.pi*np.sin(2*np.pi*t)
alpha_dot = lambda t: 0 if t<0.24 else (2*np.pi*np.deg2rad(22.5)*np.sin(2*np.pi*(t-0.25)) if t < 0.75 else 0.0)

u = lambda t: 5.0*t
h = lambda t: 0.0#1-np.cos(2*np.pi*t)
alpha = lambda t: 0 if t<0.24 else (np.deg2rad(22.5) - np.deg2rad(22.5)*np.cos(2*np.pi*(t-0.25)) if t < 0.75 else np.deg2rad(45))

be  = ao.camber_line(chords, 35, x_dot,h_dot,alpha_dot,u,h,alpha,t_step)

field   = ao.vorticity_field(chords[0])

lev_flag = 0

for t in td:

    if t > 0.0:
        print(t)
        lev_flag_prev = lev_flag

        if t > t_step:
            be.fourier_old = deepcopy(be.fourier) 

        be.update_pos(t)

        field.shed_tev(be,t)

        # be.kelvinkutta(field,0.0001,t)
        be.kelvinkutta(field,t)


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
            fig.set_size_inches(19.20, 10.80)
            plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
            contf = ax.scatter(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                    np.concatenate((field.tev_y, field.lev_y, field.ext_y))
                    , c=np.concatenate((field.tev, field.lev, field.ext)), vmin=-0.25, vmax=0.25, s=1.7)
            fig.colorbar(contf,orientation='horizontal')
            
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
            print(f'{t:.2f} {abs(be.fourier[0]) - lesp_crit:.3f} {gb}') 
        # # print(x_dot(t))

#####################################################################################    

        field.advect(be,t_step,t)

        cl[round(t/t_step)], cd[round(t/t_step)] = be.calc_cl(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                                         np.concatenate((field.tev_y, field.lev_y, field.ext_y)),
                                         np.concatenate((field.tev, field.lev, field.ext)),
                                         t, t_step)

        # print(cl[round(t/t_step)])
        # print(be.fourier[0])


fig, ax = plt.subplots()
fig.dpi = 300
ax.plot(td,
        cl,
        'k')

plt.savefig('flat1' + '.png')
plt.clf()   

np.savetxt(str(5) + 'deg.csv', np.transpose(np.vstack((cl,td))),  
          delimiter = ",")