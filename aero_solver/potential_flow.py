import math
import scipy
import scipy.integrate

'''
_N is used to denote an array of N size of a given structure, eg. Gamma_N is an array containing the vorticity strengths of N vortex blobs
'''

def eta(xi):

    return 0.0

def u_ind(Gamma_n, eta_n, xi_n, v_core, c):
    '''
    -induced velocity in the x direction due to vortex blob k
    '''

    trans = lambda theta: 0.5 * c * (1 - math.cos(theta))

    return lambda theta: 0.5 * Gamma_n / math.pi * (eta(trans(theta)) - eta_n) / math.sqrt(((trans(theta)- xi_n)**2 + (eta(trans) - eta_n)**2)**2 + v_core**4)

def v_ind(Gamma_n, eta_n, xi_n, v_core, c):
    '''
    -induced velocity in the y direction due to vortex blob k
    '''

    trans = lambda theta: 0.5 * c * (1 - math.cos(theta))

    return lambda theta: 0.5 * Gamma_n / math.pi * (trans(theta) - xi_n) / math.sqrt(((trans(theta) - xi_n)**2 + (eta(trans(theta))- eta_n)**2)**2 + v_core**4)

def dndeta(Gamma_n, eta_n, xi_n, v_core, alpha_eff, c):
    '''
    - velocity field generated by vortex blob k
    - returns dndeta(eta, xi)
    '''

    u = u_ind(Gamma_n, eta_n, xi_n, v_core, c)
    v = v_ind(Gamma_n, eta_n, xi_n, v_core, c)

    return lambda theta: u(theta) * math.sin(alpha_eff) + v(theta) * math.cos(alpha_eff)

def W_0(alpha_eff, U_ref):
    '''
    - downwash function from free stream
    - return W_0(theta, t)
    '''

    y_dot = h_dot()

    return lambda theta, t: - U_ref * math.sin(alpha_eff) + y_dot(t) * math.cos(alpha_eff)


def W_theta(W_xi, c):
    '''
    - downwash function after applying glauert transformation
    '''

    trans = lambda_xi_2_theta(c)

    return lambda theta, eta, t: W_xi(trans(theta), eta, t)

# def A_0(func,llim,ulim, eta, t):
#     '''
#     - first fourier coefficient for vorticity distribution at a given time step
#     '''
#     return 2 / math.pi * scipy.integrate.quad(func, llim, ulim, args=(eta, t))

# def A_n(func,llim,ulim, eta, t, n):
#     '''
#     - first fourier coefficient for vorticity distribution at a given time step
#     '''

#     integrand = lambda theta, eta, t: func(theta, eta, t) * math.cos(n*theta)

#     return - 1 / math.pi * scipy.integrate.quad(integrand, llim, ulim, args=(eta, t))

def Gamma_b(A_0, A_1, U_ref, c):
    '''
    - calculates the circulation of the vortex distribution on the camber line
    '''

    return math.pi * c * U_ref * (A_0 + A_1 * 0.5)

def h_dot():
    amp = 0
    p = 0
    return lambda t: amp*math.sin(p*t)

def xi_2_theta(xi,c):
    '''
    - glauert's transformation from body fixed frame
    '''

    return math.acos(1 - 2*xi / c)

def theta_2_xi(theta,c):
    '''
    - inverse of glauert's transformation
    '''

    return 0.5 * c * (1 - math.cos(theta))

def lambda_xi_2_theta(c):
    '''
    - glauert's transformation from body fixed frame
    '''

    return lambda xi: math.acos(1 - 2*xi / c)

def lambda_theta_2_xi(c):
    '''
    - inverse of glauert's transformation
    '''

    return lambda theta: 0.5 * c * (1 - math.cos(theta))

