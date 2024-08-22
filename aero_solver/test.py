import numpy as np
import matplotlib.pyplot as plt
from bem import bem


from multiprocessing import Pool

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

    return lambda t: 0.5 * np.piecewise(t, cond(t), func(t))

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

    return lambda t: 0.5 * np.piecewise(t, cond(t), func(t))

def u_func():

    return lambda t: 1 + 0.01*np.cos(0.5*np.pi*t)

# Initialise problem
U_ref = 1
alpha_eff = np.deg2rad(0)   
c = 0.5
t_step = 0.025
no_steps = 440

h1 = h_func()
h2 = h_func()

hdot1 = h_dot_func()
hdot2 = h_dot_func()

h = [h1, h2]
hdot = [hdot1, hdot2]
c = [1, 1]
alpha_eff = [0, 0]
U_ref = [1, 1]
t_step = [0.025, 0.025]
no_steps = [440, 440]

def pool_handler():
    p = Pool(2)
    p.map(bem, [U_ref, alpha_eff, c, t_step, no_steps, h, hdot])

# cl1, td1 = bem(U_ref, alpha_eff, c, t_step, no_steps, h, hdot)



pool_handler()