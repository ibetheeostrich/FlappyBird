import math
import scipy
import scipy.integrate as inte
import pot_func as pot
import pot_aux as pota
import numpy as np

import matplotlib.pyplot as plt



def main():

    # Initialise problem
    U_ref = 5
    alpha_eff = np.deg2rad(0)
    c = 4
    t_step = 0.05
    t_end = 50 * t_step
    t_d = np.arange(0,t_end,t_step)
    cl = np.array([])

    # Initialise vortex blob parameters
    N = 0
    Gamma_N = np.array([0.0])
    xi_N    = np.array([0.0])
    eta_N   = np.array([0.0])
    x_N     = np.array([0.0])
    y_N     = np.array([0.0])
    v_core  = 0.5 * c

    # Glauert's Transformation

    g_trans = lambda theta: 0.5 * c * (1 - math.cos(theta))

    # Camberline
    eta_xi = lambda xi: 0.0

    # Initialise Fourier coefficient matrix
    A = np.zeros(6)


    # Time loop

    for t in t_d:

        # print(eta_N)

        if t == t_step: # first shed vortex

            # adding vortex positions and new vortex vorticity
            xi_N[0]     = c + t*U_ref
            eta_N[0]    = 0.88*t_step
            x_N[0]      = c
            y_N[0]      = 0.0
            Gamma_N[0]  = 0.0

            N += 1


        
        # Solving for Vorticity
        Gamma_err = 10000
  
        while abs(Gamma_err) > 0.001:

            # Downwash from free stream
            W_0 = pot.W_0(U_ref, alpha_eff, t)

            # Computing Fourier coefficients
            for i in range(len(A)):

                A[i] = 0.0

                # Computing A_0 in fourier series of vorticity distribution on the bound vortex
                if i == 0: 

                    if N == 0: # solving for t=0
                        A[i], extra = inte.quad(W_0, 0.0, np.pi)
                        A[i] *= - 1.0 / np.pi / U_ref

                    else: # solving for t > 0
                            
                        A[i], extra = inte.quad(W_0, 0.0, np.pi)
                        
                        for n in range(N):                            

                            Gamma_n = Gamma_N[n]
                            eta_n   = eta_N[n]
                            xi_n    = xi_N[n]
                            dphideta = pot.dphideta(xi_n, eta_n, Gamma_n, v_core, alpha_eff)

                            integrand_n = lambda theta: dphideta(g_trans(theta))

                            A_int, extra = inte.quad(integrand_n, 0.0, np.pi)

                            A[i] += A_int

                        A[i] *= - 1.0 / np.pi / U_ref


                # Computing A_n in fourier series of vorticity distribution on the bound vortex
                else:


                    if N == 0: # solving for t=0

                        integrand_n = lambda theta: W_0(theta) * math.cos(i * theta)
                        A[i], extra = inte.quad(integrand_n , 0.0, np.pi)
                        A[i] *= 2.0 / np.pi / U_ref

                    else: # solving for t > 0

                        integrand_n = lambda theta: W_0(theta) * math.cos(i * theta)

                        A[i], extra = inte.quad(integrand_n , 0.0, np.pi)

                        for n in range(N):

                            
                            Gamma_n = Gamma_N[n]
                            eta_n   = eta_N[n]
                            xi_n    = xi_N[n]
                            dphideta = pot.dphideta(xi_n, eta_n, Gamma_n, v_core, alpha_eff)

                            integrand_n = lambda theta: dphideta(g_trans(theta)) * math.cos(i * theta)

                            A_int, extra = inte.quad(integrand_n, 0.0, np.pi)

                            A[i] += A_int

                        A[i] *= 2.0 / np.pi / U_ref
            
            Gamma_b = np.pi * c * U_ref * (A[0] + A[1] * 0.5)
            Gamma_err = Gamma_b + sum(Gamma_N)
            # print(Gamma_err, Gamma_b, sum(Gamma_N))
            Gamma_N[-1] -= 0.5*Gamma_err

        if t > 0:

            u_ind = np.zeros(N)
            v_ind = np.zeros(N)          

            for n in range(N):
                for m in range(N):

                    u_ind_p, v_ind_p = pot.V_ind_ub(xi_N[n], eta_N[n], xi_N[m], eta_N[m],  Gamma_N[m], v_core)

                    u_ind[n] += u_ind_p
                    v_ind[n] += v_ind_p

                    if m == n:
                        u_ind[n] += 0.0
                        v_ind[n] += 0.0

                trans = lambda xi: np.arccos(1 - 2*xi/c)
                gamma = lambda xi: A[0] * (1 + np.cos(trans(xi))) / np.sin(trans(xi)) + A[1] * np.sin(trans(xi)) + A[2] * np.sin(2*trans(xi)) + A[3] * np.sin(3*trans(xi))

                u_ind_p, v_ind_p = pot.V_ind_b(gamma, xi_N[n], eta_N[n], c)

                u_ind[n] = u_ind_p
                v_ind[n] = v_ind_p

            # print(u_ind, v_ind)

            xi_N     = xi_N  + u_ind*t_step #- U_ref * t_step           
            eta_N    = eta_N + v_ind*t_step - 0.88 * t_step #- 10*math.sin(1*t)* t_step  

            # print(u_ind)

            xi_N    = np.append(xi_N,(xi_N[-1] - c)*0.33 + c)
            eta_N   = np.append(eta_N,(eta_N[-1] - 0.0)*0.33 + 0)


            Gamma_N = np.append(Gamma_N, 0.0)

            N += 1


        cl = np.append(cl, - np.pi * (2 * A[0]+ A[1]))
        # print(cl, sum(Gamma_N),N)
        print(cl[-1], A[0],N)
        

    plt.plot(xi_N, 'ro')
    plt.plot([0, c],[0, 0], 'k')
    plt.axis("equal")
    plt.show()

    plt.plot(t_d,cl)
    plt.show()

    print(t_step*U_ref/c)


if __name__ == "__main__":
    main()




