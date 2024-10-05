import numpy as np
import matplotlib.pyplot as plt
from aero_solver.bem import bem

def h_dot_func():
    cond = lambda t:  [
        t < 0.5, 
        t >= 0.5 and t < 2.5,
        t >= 2.5 and t < 3.5, 
        t >= 3.5 and t < 5.5,
        t >= 5.5 and t < 6.0,
        t >= 6.0 and t < 8.0,
        t >= 8.0 and t < 9.0,
        t >= 9.0 and t < 11.0
    ]

    func = lambda t: [
        0,
        0.5*np.pi*np.sin(0.5*np.pi*(t-0.5)),
        - np.pi*2*np.sin(np.pi*(t-2.5)),
        0.5*np.pi*np.sin(0.5*np.pi*(t-3.5)),
        0,
        0.5*np.pi*np.sin(0.5*np.pi*(t-6)),
        - np.pi*2*np.sin(np.pi*(t-8)),
        0.5*np.pi*np.sin(0.5*np.pi*(t-9))
    ]

    return lambda t: np.piecewise(t, cond(t), func(t))

def h_func():
    cond = lambda t:  [
        t < 0.5, 
        t >= 0.5 and t < 2.5,
        t >= 2.5 and t < 3.5, 
        t >= 3.5 and t < 5.5,
        t >= 5.5 and t < 6.0,
        t >= 6.0 and t < 8.0,
        t >= 8.0 and t < 9.0,
        t >= 9.0 and t < 11.0
    ]

    func = lambda t: [
        0,
        1 - np.cos(0.5*np.pi*(t-0.5)),
        2*np.cos(np.pi*(t-2.5)),
        - 1 - np.cos(0.5*np.pi*(t-3.5)),
        0,
        1 - np.cos(0.5*np.pi*(t-6)),
        2*np.cos(np.pi*(t-8)),
        - 1 - np.cos(0.5*np.pi*(t-9))
    ]

    return lambda t: np.piecewise(t, cond(t), func(t))

# Initialise problem
U_ref = 5
alpha_eff = np.deg2rad(0)   
c = 2.0
t_step = 0.025
no_steps = 220

h = h_func()
hdot = h_dot_func()

cl1, td1 = bem(U_ref, alpha_eff, c, t_step, no_steps, h, hdot)

plt.plot(td1,cl1)
plt.show()
