import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from aero_solver.bem import bem as bem
from wing_kinematics.kinematics import blade_element_kinematics as bek
  
from wing_kinematics.kinematics import wing_kinematics as wk

from multiprocessing import Pool

from geom import *

import time

# Initialise problem
U_ref = 5
alpha_eff = np.deg2rad(0)   

t_step = 0.0025
no_steps = 200

no_bem = 12

# Time span
t_span = np.linspace(0.0, t_step*no_steps, no_steps, endpoint=False)


# Initialises the geometry and allows the kinematics to be solved
wing_kin = wk(I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out, F_I, F_II, A_I, A_II)

# Defining root kinematics that will drive morphing
root_kin = lambda x: - 0.05 #- 0.01 - 0.01*np.cos(2 * np.pi * 4 * x) 


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

for i in range(no_bem-1):

    le_pos_temp[i, :] = le_pos_temp[i, :] - le_pos_temp[i, 0]


for i in range(no_bem-1):

    kin = bek(np.deg2rad(40) , 4, r_pos_temp[i,:], chords_temp[i,:], le_pos_temp[i,:], U_ref,t_step)
    
    args.append((U_ref, alpha_eff, chords_temp[i,:], t_step, no_steps, kin))


# # Multiprocessing
def pool_handler():
    p = Pool(11)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    results = p.starmap(bem, args[4:5])

    return results

start = time.time()

a = pool_handler()


# print(time.time()-start)

cl_mat = np.zeros((no_bem,no_steps-1))
r_mat = r_pos_temp[:,:-2]


i = 0
for results in a:

    cl_mat[i,:] = cl = results[0]

    td = results[1]

    x = results[2]
    y = results[3]

    gamma = results[4]
    pos_arr = np.where(gamma > 0)[0]
    neg_arr = np.where(gamma < 0)[0]

    plt.plot(td,cl)
    # plt.plot(x[pos_arr],y[pos_arr],'ro')
    # plt.plot(x[neg_arr],y[neg_arr],'bo')
    # plt.axis('equal') jh    
    plt.savefig('cl v t' + str(i) + '.png',)
    plt.clf()

    i+=1

cl_int = np.trapezoid(cl_mat,r_mat,axis=0)

plt.plot(td, cl_int)
plt.savefig('clint v t' + str(i) + '.png',)
plt.clf()

# results = bem(U_ref, alpha_eff, c, t_step, no_steps, kin)

# cl = results[0]
# td = results[1]
# x = results[2]
# y = results[3]
# plt.plot(td,cl)
# plt.show()
# plt.plot(x,y,'ro')
# plt.show()


# plt.plot(td1,cl1)
# plt.show()

# plt.plot(x_N,y_N,'ro')
# plt.show()

