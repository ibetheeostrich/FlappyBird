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
U_ref = 10
alpha_eff = np.deg2rad(0)   

t_step = 0.0025
no_steps = 200

no_bem = 12

# Time span
t_span = np.linspace(0.0, t_step*no_steps, no_steps, endpoint=False)


# Initialises the geometry and allows the kinematics to be solved
wing_kin = wk(I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out, F_I, F_II, A_I, A_II)

# Defining root kinematics that will drive morphing
root_kin = lambda x: - 0.05 - 0.01 - 0.01*np.cos(2 * np.pi * 8 * x) 




# Initialising array for storing BEM kinematics

r_pos_temp = np.zeros((no_bem,no_steps+1))  
chords_temp = np.zeros((no_bem,no_steps+1))
le_pos_temp = np.zeros((no_bem,no_steps+1))
area_temp = np.zeros(no_steps+1)

for t in np.append(t_span,t_span[-1]+t_step):

    sa, w1, f1, f2, f3 = wing_kin.kin_2d_V2(t, root_kin)

    # # calculating the chord, blade element position and leading position at each time step
    # r_pos, chords, le_pos, area = wing_kin.variable_wing_params(sa, w1, f1, f2, f3, no_bem)

    # r_pos_temp[:, int(t/t_step)] = r_pos
    # chords_temp[:, int(t/t_step)] = chords
    # le_pos_temp[:, int(t/t_step)] = le_pos
    # area_temp[int(t/t_step)] = area

    x = [sa[0], w1[0], f3[0], f2[0], f1[0], root_kin(t)]
    y = [sa[1], w1[1], f3[1], f2[1], f1[1], 0.0]

    plt.plot(x,y)
    plt.axis('equal')
    plt.show()