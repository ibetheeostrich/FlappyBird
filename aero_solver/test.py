import numpy as np
import scipy.integrate as inte
import potential_flow as pot


# Initialise vortex blob parameters
N = 1
Gamma_N = [0]
eta_N = [0]
xi_N = [0]
v_core = 0

# Initialise Fourier coefficient matrix
A = [0.0, 0.0, 0.0, 0.0]

# Initialise problem
U_ref = 10.0
alpha_eff = np.deg2rad(0)
c = 1
t_step = 0.01
t_end = 10
t_d = np.arange(0,t_end,t_step)

# Testing Functions

W_0 = pot.W_0(alpha_eff, U_ref)

# Time stepping
for t in t_d:

    # Computing Fourier coefficients
    for i in range(len(A)):

        A[i] = 0.0

        # Computing A_0 in fourier series of vorticity distribution on the bound vortex
        if i == 0:

            for n in range(N):

                if n == 0:

                    integrand_0 = lambda theta, t: W_0(theta,t) / U_ref

                    A[i], extra = inte.quad(integrand_0, 0.0, np.pi,args = (t))
                    A[i]= - 1.0 / np.pi * A[i]

                else:
                    Gamma_n = Gamma_N[n]
                    eta_n   = eta_N[n]
                    xi_n    = xi_N[n]
                    dndeta = pot.dndeta(Gamma_n, eta_n, xi_n, v_core, alpha_eff, c)

                    integrand_n = lambda theta: dndeta(theta) / U_ref

                    A[i], extra = inte.quad(integrand_n, 0.0, np.pi)

                    A[i] = [i] - 1.0 / np.pi * A[i]

        # Computing A_n in fourier series of vorticity distribution on the bound vortex
        else:

            A[i] = 0.0

            for n in range(N):

                if n == 0:

                    integrand_0 = lambda theta, t: W_0(theta,t) * np.cos(theta) / U_ref

                    A[i], extra = inte.quad(integrand_0, 0.0, np.pi,args = (t))

                    A[i] = 2.0 / np.pi * A[i]

                else:
                    Gamma_n = Gamma_N[n]
                    eta_n   = eta_N[n]
                    xi_n    = xi_N[n]
                    dndeta = pot.dndeta(Gamma_n, eta_n, xi_n, v_core, alpha_eff, c)

                    integrand_n = lambda theta: dndeta(theta) * np.cos(theta) / U_ref

                    A[i], extra = inte.quad(integrand_n, 0.0, np.pi)

                    A[i] = A[i] + 2.0 / np.pi * A[i]

    # Shedding initial vortex
    if t == 0:
        eta_N[n] = U_ref * t_step + c

    # Shedding subsequent TEV
    '''
    this only does trailing edge seperation currently
    '''

    # print(A[0])
    # print(A[1])
    # print(A[2])
    # print(A[3])
    # Gamma_b = np.pi * c * U_ref * (A[0] + A[1] * 0.5)
    cl = np.pi * (2 * A[0]+ A[1])
    # print(Gamma_b)
    print(cl)

   

# Calculate circulation of bound vortex
