import math
import scipy
import scipy.integrate as inte
import pot_func as pot
import pot_aux as pota
import numpy as np
from multiprocessing import Pool
from copy import deepcopy

import matplotlib.pyplot as plt

# Constants and Globals

PI_inv = 1 / math.pi


def fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t):

    A = np.zeros(A_no)
    U_ref_inv = 1 / U_ref

    # Downwash from free stream
    W_0 = pot.W_0(U_ref, alpha_eff, t)
    # Computing Fourier coefficients
    for i in range(len(A)):
        A[i] = 0.0
        # Computing A_0 in fourier series of vorticity distribution on the bound vortex
        if i == 0: 
            if N == 0: # solving for t = 0
                A[i] = - PI_inv * U_ref_inv * pot.W_0_fast_1(U_ref, alpha_eff, t) * np.pi
            else: # solving for t > 0

                A[i] = pot.W_0_fast_1(U_ref, alpha_eff, t) * np.pi
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
            if N == 0: # solving for t = 0
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

    Gamma_b = np.pi * c * U_ref * (A[0] + A[1] * 0.5)

    return A, sum(Gamma_N), Gamma_b + sum(Gamma_N)


def main():
    # Initialise problem
    U_ref = 2
    U_ref_inv = 1/U_ref   
    alpha_eff = np.deg2rad(0)   
    c = 2.0
    t_step = 0.05
    t_end = 60 * t_step
    t_d = np.arange(0,t_end,t_step)
    cl = np.array([])


    LESP = 0.3

    LESP_flag = 0

    # Initialise vortex blob parameters
    N = 0
    Gamma_N = np.array([0.0])
    xi_N    = np.array([c])
    eta_N   = np.array([0.0])
    x_N     = np.array([c])
    y_N     = np.array([0.0])
    v_core  = 1.3*t_step*U_ref

    # Glauert's Transformation

    g_trans = lambda theta: 0.5 * c * (1 - np.cos(theta))

    # Camberline
    eta_xi = lambda xi: 0.0

    # Initialise Fourier coefficient matrix and calculation variabls
    A = np.zeros(2)
    A_no = 2

    # Newton - Raphson Params
    h = 0.01

    # Time loop

    for t in t_d:

        if t > 0:
            # Solving for TEV vorticity
            Gamma_err = 10000
  
            while abs(Gamma_err) > 0.001:

                A, Gamma_sum, Gamma_tot = fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                # print(Gamma_err, Gamma_b, sum(Gamma_N))
                Gamma_err = Gamma_tot
                Gamma_N[-1] -= Gamma_tot

                if t > 0:

                    # Newton - Raphson iteration

                    # Guess
                    x_i = deepcopy(Gamma_N[-1])

                    # inputs at guess +h and -h to estimate first derivative
                    Gamma_N_p = deepcopy(Gamma_N)
                    Gamma_N_p[-1] = x_i + h

                    Gamma_N_m = deepcopy(Gamma_N)
                    Gamma_N_m[-1] = x_i - h

                    # calculating terms for estimating first derivative
                    A, Gamma_sum, Gamma_tot_0 = fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)
                    A, Gamma_sum, Gamma_tot_p = fourier_gamma_calc(A_no, Gamma_N_p, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)
                    A, Gamma_sum, Gamma_tot_m = fourier_gamma_calc(A_no, Gamma_N_m, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                    # Newton - Raphson iteration
                    Gamma_N[-1] = x_i - Gamma_tot_0 / (0.5 * (Gamma_tot_p - Gamma_tot_m)/h)


        if abs(A[0]) > LESP:
            LESP_flag = 1
            Gamma_err = 100000

            Gamma_end = deepcopy(Gamma_N[-1])

            Gamma_N = np.append(Gamma_N, -Gamma_end)

            xi_N = np.append(xi_N, U_ref * 0.5 * t_step)
            eta_N = np.append(eta_N, - 0.5*pot.hdot(t)*t_step)

            x_N = np.append(x_N, pot.bodyin2x(xi_N[-1], t, U_ref))
            y_N = np.append(y_N, pot.bodyin2y(eta_N[-1], t))

            N += 1
                
            while abs(Gamma_err) > 0.001 and abs(abs(A[0]) - LESP) > 0.001 :

                # 2D Newton - Raphson iteration

                # Guess
                x_i = deepcopy(Gamma_N[-1])
                y_i = deepcopy(Gamma_N[-2])

                # inputs at guess +h and -h to estimate first derivative
                Gamma_N_p_LEV = deepcopy(Gamma_N)
                Gamma_N_p_LEV[-1] = x_i + h

                Gamma_N_m_LEV = deepcopy(Gamma_N)
                Gamma_N_m_LEV[-1] = x_i - h

                Gamma_N_p_TEV = deepcopy(Gamma_N)
                Gamma_N_p_TEV[-2] = y_i + h

                Gamma_N_m_TEV = deepcopy(Gamma_N)
                Gamma_N_m_TEV[-2] = y_i - h
                
                # calculating terms for estimating first derivative
                A, Gamma_sum, Gamma_tot_0 = fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)
                
                A_LEV_p, Gamma_sum, Gamma_tot_p_LEV = fourier_gamma_calc(A_no, Gamma_N_p_LEV, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)
                A_LEV_m, Gamma_sum, Gamma_tot_m_LEV = fourier_gamma_calc(A_no, Gamma_N_m_LEV, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                A_TEV_p, Gamma_sum, Gamma_tot_p_TEV = fourier_gamma_calc(A_no, Gamma_N_p_TEV, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)
                A_TEV_m, Gamma_sum, Gamma_tot_m_TEV = fourier_gamma_calc(A_no, Gamma_N_m_TEV, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                F = np.array([A[0] - LESP, Gamma_tot_0])
                J = np.array([[(A_LEV_p[0] - A_LEV_m[0]) / (2*h), (A_TEV_p[0] - A_TEV_m[0]) / (2*h)],
                              [(Gamma_tot_p_LEV - Gamma_tot_m_LEV)/(2*h), (Gamma_tot_p_TEV - Gamma_tot_m_TEV)/(2*h)]])
                
                J_inv = np.linalg.inv(J)

                [Gamma_N[-1], Gamma_N[-2]] = np.array([x_i, y_i]) - J_inv@F 

        # Advecting and shedding vortices for next time step
        if t == 0:
            xi_N    = np.array([c + U_ref*t_step, c + U_ref*t_step/3])
            eta_N   = np.array([pot.h(t_step), pot.h(t_step)/3])

            x_N = np.array([c , c - U_ref*t_step/3])
            y_N = np.array([0, -pot.h(t_step)/3])

            Gamma_N = np.append(Gamma_N, 0.0)

            N += 2

        if t > 0:


            u_ind = np.zeros(N)
            v_ind = np.zeros(N)    


            for n in range(len(u_ind)):
                for m in range(len(u_ind)):

                    u_ind_p, v_ind_p = pot.V_ind_ub(x_N[n], y_N[n], x_N[m], y_N[m],  Gamma_N[m], v_core)

                    u_ind[n] += u_ind_p
                    v_ind[n] += v_ind_p

                    if m == n:
                        u_ind[n] += 0.0
                        v_ind[n] += 0.0
#################################################################################################################################################################################################
#              MAYBE NEED?????????               
#################################################################################################################################################################################################                
                trans = lambda xi: np.arccos(1 - 2*xi/c)
                gamma = lambda xi: A[0] * (1 + np.cos(trans(xi))) / np.sin(trans(xi)) + A[1] * np.sin(trans(xi)) #+ A[2] * np.sin(2*trans(xi)) + A[3] * np.sin(3*trans(xi))

                u_ind_p, v_ind_p = pot.V_ind_b(gamma, xi_N[n], eta_N[n], c)

                u_ind[n] += u_ind_p # - U_ref
                v_ind[n] += v_ind_p # - pot.hdot(t) 

            # print(u_ind, v_ind)
#################################################################################################################################################################################################
 

            x_N     = x_N + u_ind*t_step 
            y_N     = y_N + v_ind*t_step 


            if t == t_step:
                x_N = np.append(x_N,(x_N[0] - (c-U_ref*t_step))*0.33 + c-U_ref*t_step)
                y_N = np.append(y_N,(y_N[0] - pot.h(t))*0.33 + pot.h(t))
            else:
                if LESP_flag:
                    x_N = np.append(x_N,(x_N[-2] - (c-U_ref*t))*0.33 + c-U_ref*t)
                    y_N = np.append(y_N,(y_N[-2] - pot.h(t))*0.33 + pot.h(t))
                else:
                    x_N = np.append(x_N,(x_N[-1] - (c-U_ref*t))*0.33 + c-U_ref*t)
                    y_N = np.append(y_N,(y_N[-1] - pot.h(t))*0.33 + pot.h(t))



            xi_N    = pot.xin2body(x_N, t, U_ref)
            eta_N   = pot.yin2body(y_N, t)


            Gamma_N = np.append(Gamma_N, 0.0)

            N += 1

            LESP_flag = 0


        cl = np.append(cl, np.pi * (2 * A[0]+ A[1]))
        # print(cl, sum(Gamma_N),N)
        V = pot.hdot(t)
        # print(cl[-1], np.rad2deg(np.arctan2(V,U_ref)), Gamma_N[-1],A[0],t)
        print(N,t, A[0])
        # plt.plot(x_N, y_N, 'ro')
        # plt.plot([0.0-U_ref *(t), c-U_ref*(t)], [pot.h((t)), pot.h((t))], 'k')
        # plt.axis("equal")
        # plt.show()

    print(len(Gamma_N))
    plt.plot(x_N, y_N, 'ro')
    plt.plot([0.0-U_ref *(t_end-t_step), c-U_ref*(t_end-t_step)], [pot.h(t_end-t_step), pot.h(t_end-t_step)], 'k')
    plt.axis("equal")
    plt.show()

    # plt.plot(t_d,cl)
    # plt.show()
    


    

    print(t_step*U_ref/c)


if __name__ == "__main__":
    main()




