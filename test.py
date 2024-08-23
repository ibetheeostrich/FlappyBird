import numpy as np
import matplotlib.pyplot as plt
from aero_solver.bem import bem as bem
from wing_kinematics.kinematics import blade_element_kinematics as bek

from wing_kinematics.kinematics import wing_kinematics as wk

from multiprocessing import Pool

from geom import *

# Initialise problem
U_ref = 1
alpha_eff = np.deg2rad(0)   
c = 0.5
t_step = 0.025
no_steps = 220

chord = np.zeros(441) + c
le_pos = np.zeros(441)

freq = 0.4

amplitude = 0.01
span_pos = 1

kin = bek(amplitude, freq, span_pos, chord, le_pos, U_ref)

# Initialising blade element cases
args = (U_ref, alpha_eff, c, t_step, no_steps, kin)


wing_kin = wk(I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out, F_I, F_II, A_I, A_II)

root_kin = lambda x: - 0.05 - 0.01 + 0.01*np.cos(x)

ttest = np.linspace(0.0, 2*np.pi, 10)

for t in ttest:

    sa, e2, w1, f1, f2, f3 = wing_kin.kin_2d(t, root_kin)

    plt.plot([sa[0], e2[0], w1[0]], [sa[1], e2[1], w1[1]],'b')
    plt.plot([w1[0], f1[0]],[w1[1], f1[1]], 'r')
    plt.plot([w1[0], f2[0]],[w1[1], f2[1]], 'r')
    plt.plot([w1[0], f3[0]],[w1[1], f3[1]], 'r')


    plt.axis('equal')

    plt.show()


# Multiprocessing
def pool_handler():
    p = Pool(12)
    results = p.starmap(bem, [args, args, args, args, args, args, args, args, args, args, args, args])

    return results




# a = pool_handler()



# plt.plot(td1,cl1)
# plt.show()

# plt.plot(x_N,y_N,'ro')
# plt.show()
 