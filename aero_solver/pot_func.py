import math
import scipy.integrate as inte
import numpy as np
import numpy as np

PI = math.pi
PI_inv = 1 / math.pi


def hdot(t):
    if t < 0.5:
        return 0
    elif t < 2.5:
        return 0.5*np.pi*np.sin(0.5*np.pi*(t-0.5))
    elif t < 3.5:
        return - np.pi*2*np.sin(np.pi*(t-2.5))
    elif t < 5.5:
        return 0.5*np.pi*np.sin(0.5*np.pi*(t-3.5))
    
    elif t < 6:
        return 0
    elif t < 8:
        return 0.5*np.pi*np.sin(0.5*np.pi*(t-6))
    elif t < 9:
        return - np.pi*2*np.sin(np.pi*(t-8))
    elif t < 11:
        return 0.5*np.pi*np.sin(0.5*np.pi*(t-9))

def h(t):
    if t < 0.5:
        return 0
    elif t < 2.5:
        return 1 - np.cos(0.5*np.pi*(t-0.5))
    elif t < 3.5:
        return 2*np.cos(np.pi*(t-2.5))
    elif t < 5.5:
        return - 1 - np.cos(0.5*np.pi*(t-3.5))
    
    elif t < 6:
        return 0
    elif t < 8:
        return 1 - np.cos(0.5*np.pi*(t-6))
    elif t < 9:
        return 2*np.cos(np.pi*(t-8))
    elif t < 11:
        return - 1 - np.cos(0.5*np.pi*(t-9))

# def hdot(t):
#     return 1.0

# def h(t):
#     return t

# def hdot(t):
#     return 0

# def h(t):
#     return 0



def xin2body(x,t, U_ref):
    return x + U_ref*t

def yin2body(y,t):
    return y - h(t)

def bodyin2x(x,t, U_ref):
    return x - U_ref*t

def bodyin2y(y,t):
    return y + h(t)


def g_trans(theta, c):
    return 0.5 * c * (1 - np.cos(theta))

def dphideta(xi_n, eta_n , Gamma_n, v_core, alpha_eff):

    const = 0.5 * Gamma_n * PI_inv

    eta = lambda xi: 0.0

    func = lambda xi:  const * (
        (-eta_n * math.sin(alpha_eff) - (xi - xi_n) * math.cos(alpha_eff))  
        / 
        np.sqrt((eta_n**2 + (xi - xi_n)**2)**2 + v_core**4)
                         )
    return func

def W_0(U_ref, alpha_eff, t):

    A = 1
    p = 0.5

    eta = lambda xi: 0.0

    return lambda xi: - U_ref*math.sin(alpha_eff) + (hdot(t)) * math.cos(alpha_eff)

def W_0_fast_1(U_ref, alpha_eff, t):

    return - U_ref*math.sin(alpha_eff) + (hdot(t)) * math.cos(alpha_eff)

def V_ind_b(gamma, xi_n, eta_n, c):
    '''
    - takes in vorticity distribution gamma
    - find the induced velocity of the vorticity distribution at point (xi_n, eta_n)
    '''

    integrand_u = lambda xi: gamma(xi) * (0.0 - eta_n) / ((xi_n - xi)**2 + (eta_n - 0.0)**2)

    def_int_u, extra = inte.quad(integrand_u, 0, c)

    u_ind = 0.5 * PI_inv * def_int_u

    integrand_v = lambda xi: gamma(xi) * (xi - xi_n) / ((xi_n - xi)**2 + (eta_n - 0.0)**2)

    def_int_v, extra = inte.quad(integrand_v, 0, c)

    v_ind = -0.5 * PI_inv * def_int_v 

    return u_ind, v_ind

def V_ind_ub(xi, eta, xi_n, eta_n, gamma, v_core):

    '''
    - calculates the induced velocity at a point by another vortex blob
    '''

    u_ind_ub = 0.5 * gamma * PI_inv * (eta - eta_n) / math.sqrt(((xi - xi_n)**2 + (eta - eta_n)**2)**2 + v_core**4)

    v_ind_ub = -0.5 * gamma * PI_inv * (xi - xi_n) / math.sqrt(((xi - xi_n)**2 + (eta - eta_n)**2)**2 + v_core**4)

    return u_ind_ub, v_ind_ub

