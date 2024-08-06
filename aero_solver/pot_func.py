import math

PI = math.pi
PI_inv = 1 / math.pi




def dphideta(xi_n, eta_n , Gamma_n, v_core, alpha_eff):

    const = 0.5 * Gamma_n * PI_inv

    eta = lambda xi: 0.0

    func = lambda xi:  const * (
        ((eta(xi) - eta_n) * math.sin(alpha_eff) + (xi - xi_n) * math.cos(alpha_eff))  
        / 
        (((eta(xi) - eta_n)**2 + (xi - xi_n)**2)**2 + v_core**4)
                         )


    return func

def W_0(U_ref, alpha_eff, t):

    A = 1
    p = 0.5

    eta = lambda xi: 0.0

    return lambda xi: - U_ref*math.sin(alpha_eff) + (eta(xi) - A*math.sin(p* t)) * math.cos(alpha_eff)