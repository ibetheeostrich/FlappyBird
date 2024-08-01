import numpy as np
import scipy.integrate as inte
import potential_flow as pot


# Initialise vortex blob parameters
N = 1
Gamma_N = [0]
eta_N = [0]
xi_N = [0]
v_core = 0

# Initialise problem
U_ref = 10.0
alpha_eff = np.deg2rad(10)
c = 1

# Testing Functions

W_0 = pot.W_0(alpha_eff, U_ref)


# Computing A_0 in fourier series of vorticity distribution on the bound vortex
A_0 = 0

for n in range(N):

    if n == 0:

        integrand_0 = lambda theta, t: W_0(theta,t  ) / U_ref

        A_0_int = inte.quad(integrand_0, 0.0, np.pi,args = (0.0))
        A_0 = - 1.0 / np.pi * A_0_int[0]

    else:
        Gamma_n = Gamma_N[n]
        eta_n   = eta_N[n]
        xi_n    = xi_N[n]
        dndeta = pot.dndeta(Gamma_n, eta_n, xi_n, v_core, alpha_eff, c)

        integrand_n = lambda theta: dndeta(theta) / U_ref

        A_0_int = inte.quad(integrand_n, 0.0, np.pi)

        A_0 = A_0 - 1.0 / np.pi * A_0_int[0]

# Computing A_1 in fourier series of vorticity distribution on the bound vortex
A_1 = 0

for n in range(N):

    if n == 0:

        integrand_0 = lambda theta, t: W_0(theta,t) * np.cos(theta) / U_ref

        A_1_int = inte.quad(integrand_0, 0.0, np.pi,args = (0.0))

        A_1 = 2.0 / np.pi * A_1_int[0]

    else:
        Gamma_n = Gamma_N[n]
        eta_n   = eta_N[n]
        xi_n    = xi_N[n]
        dndeta = pot.dndeta(Gamma_n, eta_n, xi_n, v_core, alpha_eff, c)

        integrand_n = lambda theta: dndeta(theta) * np.cos(theta) / U_ref

        A_1_int = inte.quad(integrand_n, 0.0, np.pi)

        A_1 = A_1 + 2.0 / np.pi * A_1_int[0]

# Calculate circulation of bound vortex
Gamma_b = np.pi * c * U_ref * (A_0 + A_1 * 0.5)

cl = np.pi * (2 * A_0+ A_1)

print(A_0)
print(A_1)
print(Gamma_b)
print(cl)