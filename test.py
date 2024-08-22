import numpy as np
import matplotlib.pyplot as plt
from aero_solver.bem import bem as bem
from wing_kinematics.kinematics import blade_element_kinematics as bek

from multiprocessing import Pool

# def h_dot_func():
#     cond = lambda t:  [
#         t < 0.5, 
#         t >= 0.5 and t < 2.5,
#         t >= 2.5 and t < 3.5, 
#         t >= 3.5 and t < 5.5,
#         t >= 5.5 and t < 6.0,
#         t >= 6.0 and t < 8.0,
#         t >= 8.0 and t < 9.0,
#         t >= 9.0 and t < 11.0
#     ]

#     func = lambda t: [
#         0,
#         0.5*np.pi*np.sin(0.5*np.pi*(t-0.5)),
#         - np.pi*2*np.sin(np.pi*(t-2.5)),
#         0.5*np.pi*np.sin(0.5*np.pi*(t-3.5)),
#         0,
#         0.5*np.pi*np.sin(0.5*np.pi*(t-6)),
#         - np.pi*2*np.sin(np.pi*(t-8)),
#         0.5*np.pi*np.sin(0.5*np.pi*(t-9))
#     ]

#     return lambda t: 0.5 * np.piecewise(t, cond(t), func(t))

# def h_func():
#     cond = lambda t:  [
#         t < 0.5, 
#         t >= 0.5 and t < 2.5,
#         t >= 2.5 and t < 3.5, 
#         t >= 3.5 and t < 5.5,
#         t >= 5.5 and t < 6.0,
#         t >= 6.0 and t < 8.0,
#         t >= 8.0 and t < 9.0,
#         t >= 9.0 and t < 11.0
#     ]

#     func = lambda t: [
#         0,
#         1 - np.cos(0.5*np.pi*(t-0.5)),
#         2*np.cos(np.pi*(t-2.5)),
#         - 1 - np.cos(0.5*np.pi*(t-3.5)),
#         0,
#         1 - np.cos(0.5*np.pi*(t-6)),
#         2*np.cos(np.pi*(t-8)),
#         - 1 - np.cos(0.5*np.pi*(t-9))
#     ]

#     return lambda t: 0.5 * np.piecewise(t, cond(t), func(t))

# def u_func():

#     return lambda t: 1 + 0.01*np.cos(0.5*np.pi*t)


# h1 = h_func()
# h2 = h_func()

# hdot1 = h_dot_func()
# hdot2 = h_dot_func()

# h = [h1, h2]
# hdot = [hdot1, hdot2]
# c = [1, 1]
# alpha_eff = [0, 0]
# U_ref = [1, 1]
# t_step = [0.025, 0.025]
# no_steps = [440, 440]

# def pool_handler():
#     p = Pool(2)
#     p.map(bem, [U_ref, alpha_eff, c, t_step, no_steps, h, hdot])



# Initialise problem
U_ref = 1
alpha_eff = np.deg2rad(0)   
c = 0.5
t_step = 0.025
no_steps = 220

chord = np.zeros(441) + c
le_pos = np.zeros(441)

freq = 0.4

amplitude = 1
span_pos = 1

kin = bek(amplitude, freq, span_pos, chord, le_pos, U_ref)

cl1, td1, x_N, y_N = bem(U_ref, alpha_eff, c, t_step, no_steps, kin)



t_span = np.linspace(0,t_step*no_steps,no_steps)
hdot = np.zeros(no_steps)

# i = 0
# for t in t_span:
#     hdot[i] =  kin.h_dot(t)
#     i+=1

# plt.plot(t_span,hdot)

plt.plot(td1,cl1)


plt.show()

plt.plot(x_N,y_N,'ro')
plt.show()
