import numpy as np
import scipy.integrate as inte
import potential_flow as pot
import copy


# Initialise vortex blob parameters
N = 0
Gamma_N = [0]
xi_N = [0]
eta_N = [0]
x_N = [0]
y_N = [0]
v_core = 0.001

# Camberline
eta_xi = lambda xi: 0.0

# Initialise Fourier coefficient matrix
A = [0.0, 0.0, 0.0, 0.0]

# Initialise problem
U_ref = 10.0
alpha_eff = np.deg2rad(0)
c = 1
t_step = 0.01
t_end = 100*t_step
t_d = np.arange(0,t_end,t_step)

# Testing Functions

W_0 = pot.W_0(alpha_eff, U_ref)

# Time stepping
for t in t_d:


    # Advecting shed vortices

    # Shedding initial vortex
    if t == 0:
        N = 1
        xi_N[0] = c
        eta_N[0] = 0
        
        x_N[0] = c
        y_N[0] = 1

    # Advecting first zero vortex and shedding second vortex
    if t == t_step:

        xi_N[0]     = pot.trans_x2xi(x_N[0], U_ref, t)
        eta_N[0]    = pot.trans_y2eta(y_N[0], t)

        # generating new vortex
        Gamma_N.append(1.0)


        xi_N.append((xi_N[0] - c)*0.33 + c)
        eta_N.append((eta_N[0] - 0.0)*0.33 + 0)

        x_N.append(pot.trans_xi2x(xi_N[-1], t, U_ref))
        y_N.append(pot.trans_eta2y(eta_N[-1], t))

        N += 1

    # Advecting vortices after the initial cases
    elif t > t_step:



        u_ind = np.zeros(N)
        v_ind = np.zeros(N)
        for n in range(N):

            for m in range(N):

                u_ind_p, v_ind_p = pot.V_ind_ub(xi_N[n], eta_N[n], xi_N[m], eta_N[m],  Gamma_N[m], v_core)

                u_ind[n] += u_ind_p
                v_ind[n] += v_ind_p


            trans = lambda xi: np.arccos(1 - 2*xi/c)
            gamma = lambda xi: A[0] * (1 + np.cos(trans(xi))) / np.sin(trans(xi)) + A[1] * np.sin(trans(xi)) + A[2] * np.sin(2*trans(xi)) + A[3] * np.sin(trans(xi))

            u_ind_p, v_ind_p = pot.V_ind_b(gamma, xi_N[n], eta_N[n], c, eta_xi)
            
            u_ind[n] += u_ind_p
            v_ind[n] += v_ind_p

        xi_N += u_ind * t_step
        eta_N += v_ind * t_step

        Gamma_N.append(1.0)

        xi_N = np.append(xi_N, (xi_N[0] - c)*0.33 + c)
        eta_N = np.append(eta_N, (eta_N[0] - 0.0)*0.33 + 0)  + 2 * np.sin(2*t) * t_step

        N += 1   

    
    # Computing vortex strenghts

    if t > 0:
        Gamma_err = 10000

        while abs(Gamma_err) > 0.000001:

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

                            A[i] = A[i] - 1.0 / np.pi * A[i]

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


            Gamma_b = np.pi * c * U_ref * (A[0] + A[1] * 0.5)
            Gamma_err = Gamma_b - sum(Gamma_N)
            Gamma_N[-1] += 0.5*Gamma_err

            # print(Gamma_err)


#     cl = np.pi * (2 * A[0]+ A[1])
#     print(cl)



# print(sum(Gamma_N))

    # print(A[0])
    # print(A[1])
    # print(A[2])
    # print(A[3])

    # cl = np.pi * (2 * A[0]+ A[1])
    # # print(Gamma_b)
    # print(cl)


