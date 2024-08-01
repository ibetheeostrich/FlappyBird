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

        A_1, extra = inte.quad(integrand_0, 0.0, np.pi,args = (0.0))

        A_1 = 2.0 / np.pi * A_1

    else:
        Gamma_n = Gamma_N[n]
        eta_n   = eta_N[n]
        xi_n    = xi_N[n]
        dndeta = pot.dndeta(Gamma_n, eta_n, xi_n, v_core, alpha_eff, c)

        integrand_n = lambda theta: dndeta(theta) * np.cos(theta) / U_ref

        A_1, extra = inte.quad(integrand_n, 0.0, np.pi)

        A_1 = A_1 + 2.0 / np.pi * A_1

# Computing A_1 in fourier series of vorticity distribution on the bound vortex
A_2 = 0

for n in range(N):

    if n == 0:

        integrand_0 = lambda theta, t: W_0(theta,t) * np.cos(2*theta) / U_ref

        A_2, extra = inte.quad(integrand_0, 0.0, np.pi,args = (0.0))

        A_2 = 2.0 / np.pi * A_2

    else:
        Gamma_n = Gamma_N[n]
        eta_n   = eta_N[n]
        xi_n    = xi_N[n]
        dndeta = pot.dndeta(Gamma_n, eta_n, xi_n, v_core, alpha_eff, c)

        integrand_n = lambda theta: dndeta(theta) * np.cos(2*theta) / U_ref

        A_2, extra = inte.quad(integrand_n, 0.0, np.pi)

        A_2 = A_2 + 2.0 / np.pi * A_2

# Computing A_1 in fourier series of vorticity distribution on the bound vortex
A_3 = 0

for n in range(N):

    if n == 0:

        integrand_0 = lambda theta, t: W_0(theta,t) * np.cos(3*theta) / U_ref

        A_3, extra = inte.quad(integrand_0, 0.0, np.pi,args = (0.0))

        A_3 = 2.0 / np.pi * A_3

    else:
        Gamma_n = Gamma_N[n]
        eta_n   = eta_N[n]
        xi_n    = xi_N[n]
        dndeta = pot.dndeta(Gamma_n, eta_n, xi_n, v_core, alpha_eff, c)

        integrand_n = lambda theta: dndeta(theta) * np.cos(3*theta) / U_ref

        A_3, extra = inte.quad(integrand_n, 0.0, np.pi)

        A_3 = A_3 + 2.0 / np.pi * A_3

# Calculate circulation of bound vortex
Gamma_b = np.pi * c * U_ref * (A_0 + A_1 * 0.5)

cl = np.pi * (2 * A_0+ A_1)

print(A_0)
print(A_1)
print(A_2)
print(A_3)

print(Gamma_b)
print(cl)