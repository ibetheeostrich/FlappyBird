import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from aero_solver.bem import bem as bem
from wing_kinematics.kinematics import blade_element_kinematics as bek
  
from wing_kinematics.kinematics import wing_kinematics as wk

from multiprocessing import Pool

from geom import *

import time

import os

start = time.time()


# Initialise problem
rho = 1.225
U_ref = 8
alpha_eff = np.deg2rad(4)   
wing_amplitude = 42.5

# t_step = 0.0025
no_steps = 200

no_bem = 12 

frequency = 4

t_end = 3/frequency
t_step = t_end/no_steps


# Time span
# t_span = np.linspace(0.0, t_step*no_steps, no_steps, endpoint=False)
t_span = np.linspace(0.0, t_end, no_steps, endpoint=False)


# Initialises the geometry and allows the kinematics to be solved
wing_kin = wk(I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out, F_I, F_II, A_I, A_II)

# Defining root kinematics that will drive morphing
root_kin = lambda x: - 0.05 - 0.01 - 0.01*np.cos(2 * np.pi * frequency * x) 

       
# Initialising array for storing BEM kinematics

r_pos_temp = np.zeros((no_bem,no_steps+1))  
chords_temp = np.zeros((no_bem,no_steps+1))
le_pos_temp = np.zeros((no_bem,no_steps+1))
area_temp = np.zeros(no_steps+1)

for t in np.append(t_span,t_span[-1]+t_step):

    sa, w1, f1, f2, f3 = wing_kin.kin_2d_V2(t, root_kin)


    # calculating the chord, blade element position and leading position at each time step
    r_pos, chords, le_pos, area = wing_kin.variable_wing_params(sa, w1, f1, f2, f3, no_bem)

    r_pos_temp[:, round(t/t_step)] = r_pos
    chords_temp[:, round(t/t_step)] = chords
    area_temp[round(t/t_step)] = area
    le_pos_temp[:, round(t/t_step)] = le_pos

args = []

tag = np.linspace(0,no_bem,no_bem,endpoint=False)

for i in range(no_bem-1):

    le_pos_temp[i, :] = le_pos_temp[i, :] - le_pos_temp[i, 0]


for i in range(no_bem-1):

    kin = bek(np.deg2rad(wing_amplitude) , frequency, r_pos_temp[i,:], chords_temp[i,:], le_pos_temp[i,:], U_ref, alpha_eff, t_step)
    
    args.append((round(tag[i]),U_ref, alpha_eff, chords_temp[i,:], t_step, no_steps, kin))


# # Multiprocessing
def pool_handler():
    p = Pool(no_bem-1)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    results = p.starmap(bem, args)

    return results

# Evaluating blade elements using multi processing
a = pool_handler()

# Results accumulation matrices 
cl_mat = np.zeros((no_bem,no_steps-1))
r_mat = r_pos_temp[:,1:-1]

# Processing results for each blade element
for results in a:

    # Getting BEM tag
    i = results[0]

    # Formatting Data
    cl_mat[i,:] = cl = results[1]

    td = results[2]

    x = results[3]
    y = results[4]

    gamma = results[5]
    pos_arr = np.where(gamma > 0)[0]
    neg_arr = np.where(gamma < 0)[0]

    # # Creating Figures
    # fig, ax = plt.subplots()
    # fig.dpi = 300
    # ax.plot(td, cl)
    # ax.set_xlabel('Time  (s)')
    # ax.set_ylabel('Lift Force (n)')
    # plt.savefig('l v t' + str(i) + '.png')
    # plt.close(fig)

# Integrating BEM
l_int = np.trapz(cl_mat,r_mat,axis=0)

print(np.trapz(l_int,td)) 

fig, ax = plt.subplots()
fig.dpi = 300
ax.plot(td, l_int)
ax.set_xlabel('Time  (s)')
ax.set_ylabel('Lift Force (n)')
plt.savefig('lint v t' + str(i) + '.png')
plt.show()

print(time.time()-start)