import math

def cart_vortex(x,y):

    u = 0.5 / math.pi() / (x**2 + y**2) * y 
    v = 0.5 / math.pi() / (x**2 + y**2) * x

    return u, v

def phi(xi,eta):
    return

def dphideta(xi,eta,n):

    for i in range(n):
        func = 0

    return

def h_dot(t):
    return

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

def W_1(xi, eta, xi_n, eta_n, Gamma_n):
    '''
    - is the downwash function at each position on the camber line due to shed vortices
    - used for setting the boundary condition on the camber line
    '''
    return (0.5*Gamma_n/math.pi()) * (xi - xi_n) / (((eta-eta_n)**2)*(((xi-xi_n)**2)/((eta-eta_n)**2)+  1))
    
def W_2(alpha_eff, U_ref):
    '''
    - "down wash" on the wing due to the free stream
    '''
    return U_ref*math.sin(alpha_eff)

def W_3(eta,t,alpha_eff, A):
    '''
    - "down wash" due to the vertical motion of the wing
    '''
    return (eta - A*math.sin(t))*math.cos(alpha_eff)