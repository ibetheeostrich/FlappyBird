import numpy as np
import matplotlib.pyplot as plt
from aero_solver.bem import bem as bem
from wing_kinematics.kinematics import blade_element_kinematics as bek

from wing_kinematics.kinematics import wing_kinematics as wk

from multiprocessing import Pool

from geom import *

import time

# Initialise problem
U_ref = 1
alpha_eff = np.deg2rad(0)   
c = 0.5
t_step = 0.025
no_steps = 100

no_bem = 12

# Time span
t_span = np.linspace(0.0, t_step*no_steps, no_steps)


# Initialises the geometry and allows the kinematics to be solved
wing_kin = wk(I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out, F_I, F_II, A_I, A_II)

# Defining root kinematics that will drive morphing
root_kin = lambda x: - 0.05 - 0.01 - 0.01*np.cos(2 * np.pi * 1 * x) if x >= 1 else -0.07




# Initialising array for storing BEM kinematics

r_pos_temp = np.zeros((no_bem,no_steps+1))
chords_temp = np.zeros((no_bem,no_steps+1))
le_pos_temp = np.zeros((no_bem,no_steps+1))
area_temp = np.zeros(no_steps+1)

for t in np.append(t_span,t_span[-1]+t_step):

    sa, w1, f1, f2, f3 = wing_kin.kin_2d(t, root_kin)

    # calculating the chord, blade element position and leading position at each time step
    r_pos, chords, le_pos, area = wing_kin.variable_wing_params(sa, w1, f1, f2, f3, no_bem)

    r_pos_temp[:, int(t/t_step)] = r_pos
    chords_temp[:, int(t/t_step)] = chords
    le_pos_temp[:, int(t/t_step)] = le_pos
    area_temp[int(t/t_step)] = area


args = []

for i in range(no_bem):

    kin = bek(np.deg2rad(30) , 1, r_pos_temp[i,:], chords_temp[i,:], le_pos_temp[i,:], U_ref,t_step)
    
    args.append((U_ref, alpha_eff, c, t_step, no_steps, kin))

# # Multiprocessing
def pool_handler():
    p = Pool(12)
    results = p.starmap(bem, args)

    return results




a = pool_handler()

for results in a:

    cl = results[0]
    td = results[1]

    x = results[2]
    y = results[3]

    plt.plot(td,cl)
    plt.show()

    plt.plot(x,y,'ro')
    plt.show()



print(a)



# plt.plot(td1,cl1)
# plt.show()

# plt.plot(x_N,y_N,'ro')
# plt.show()
 