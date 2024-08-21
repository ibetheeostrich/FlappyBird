import math
import scipy
import scipy.integrate as inte
import pot_func as pot
import pot_aux as pota
import numpy as np
from multiprocessing import Pool
from copy import deepcopy
import time

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import io

# Constants and Globals

PI_inv = 1 / math.pi


def fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t):

    x = np.linspace(0.0, np.pi, 513, endpoint=True)

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
                A[i] = - PI_inv * U_ref_inv * pot.W_0_fast_1(U_ref, alpha_eff, t) * c
            else: # solving for t > 0

                A[i] = pot.W_0_fast_1(U_ref, alpha_eff, t) * c

                for n in range(N):                            
                    Gamma_n = Gamma_N[n]
                    eta_n   = eta_N[n]
                    xi_n    = xi_N[n]
                    dphideta = pot.dphideta(xi_n, eta_n, Gamma_n, v_core, alpha_eff)
                    integrand_n = lambda theta: dphideta(g_trans(theta))

                    A_int = inte.trapz(integrand_n(x), x)

                    A[i] -= A_int
                A[i] *= - 1.0 / np.pi / U_ref
        # Computing A_n in fourier series of vorticity distribution on the bound vortex
        else:
            if N == 0: # solving for t = 0
                integrand_n = lambda theta: pot.W_0_fast_1(U_ref, alpha_eff, t) * math.cos(i * theta)
                A[i], extra = inte.quad(integrand_n , 0.0, np.pi)
                A[i] *= 2.0 / np.pi / U_ref
            else: # solving for t > 0
                integrand_n = lambda theta: pot.W_0_fast_1(U_ref, alpha_eff, t) * np.cos(i * theta)
                A[i], extra = inte.quad(integrand_n , 0.0, np.pi)
                for n in range(N):

                    Gamma_n = Gamma_N[n]
                    eta_n   = eta_N[n]
                    xi_n    = xi_N[n]
                    dphideta = pot.dphideta(xi_n, eta_n, Gamma_n, v_core, alpha_eff)
                    integrand_n = lambda theta: dphideta(g_trans(theta)) * np.cos(i * theta)

                    A_int = inte.trapz(integrand_n(x), x)

                    A[i] -= A_int

                A[i] *= 2.0 / np.pi / U_ref

    Gamma_b = np.pi * c * U_ref * (A[0] + A[1] * 0.5)

    return A, sum(Gamma_N), Gamma_b + sum(Gamma_N)


def main():
    start = time.time()
    # movie
    frames = []

    # Initialise problem
    U_ref = 5
    U_ref_inv = 1/U_ref   
    alpha_eff = np.deg2rad(0)   
    c = 2.0
    t_step = 0.025
    t_end = 440 * t_step
    t_d = np.arange(0,t_end,t_step)
    cl = np.array([])


    LESP = 0.2

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
    A_no = 4

    # Newton - Raphson Params
    h = 0.01

    # Time loop

    for t in t_d:

        # TEV Shedding
        if t > 0:
            # Solving for TEV vorticity
            Gamma_err = 10000
  
            while abs(Gamma_err) > 0.00001:

                A, Gamma_sum, Gamma_tot = fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                # print(Gamma_err, Gamma_b, sum(Gamma_N))
                Gamma_err = Gamma_tot

                if t > 0:

                    # Newton - Raphson iteration

                    # Guess
                    x_i = deepcopy(Gamma_N[-1])

                    # inputs at guess +h and -h to estimate first derivative
                    # calculating terms for estimating first derivative
                    Gamma_N_p = Gamma_N_m = Gamma_N

                    Gamma_N_p[-1] = x_i + h
                    A, Gamma_sum, Gamma_tot_p = fourier_gamma_calc(A_no, Gamma_N_p, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                    Gamma_N_m[-1] = x_i - h
                    A, Gamma_sum, Gamma_tot_m = fourier_gamma_calc(A_no, Gamma_N_m, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)


                    # b = Gamma_tot - Gamma_sum

                    # print(b)

                    # Newton - Raphson iteration
                    Gamma_N[-1] = x_i - Gamma_tot / (0.5 * (Gamma_tot_p - Gamma_tot_m)/h)

        # LEV Shedding
        if abs(A[0]) > LESP:
            LESP_flag = 1
            Gamma_err = 100000

            Gamma_end = deepcopy(Gamma_N[-1])

            Gamma_N = np.append(Gamma_N, -Gamma_end)

            xi_N = np.append(xi_N, 0)
            eta_N = np.append(eta_N, 0)

            x_N = np.append(x_N, pot.bodyin2x(xi_N[-1], t, U_ref))
            y_N = np.append(y_N, pot.bodyin2y(eta_N[-1], t))

            N += 1
                
            while abs(Gamma_err) > 0.001 and abs(abs(A[0]) - LESP) > 0.001 :

                A, Gamma_sum, Gamma_tot_0 = fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                # 2D Newton - Raphson iteration
                # Guess
                x_i = deepcopy(Gamma_N[-1])
                y_i = deepcopy(Gamma_N[-2])

                # inputs at guess +h and -h to estimate first derivative
                Gamma_N_p_LEV = Gamma_N_m_LEV = Gamma_N_p_TEV = Gamma_N_m_TEV =      Gamma_N

                Gamma_N_p_LEV[-1] = x_i + 4*h
                Gamma_N_p_LEV[-2] = y_i
                A_LEV_p, Gamma_sum, Gamma_tot_p_LEV = fourier_gamma_calc(A_no, Gamma_N_p_LEV, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                Gamma_N_m_LEV[-1] = x_i - 4*h
                Gamma_N_m_LEV[-2] = y_i
                A_LEV_m, Gamma_sum, Gamma_tot_m_LEV = fourier_gamma_calc(A_no, Gamma_N_m_LEV, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                Gamma_N_p_TEV[-2] = y_i + 4*h
                Gamma_N_p_TEV[-1] = x_i
                A_TEV_p, Gamma_sum, Gamma_tot_p_TEV = fourier_gamma_calc(A_no, Gamma_N_p_TEV, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

                Gamma_N_m_TEV[-2] = y_i - 4*h
                Gamma_N_m_TEV[-1] = x_i
                A_TEV_m, Gamma_sum, Gamma_tot_m_TEV = fourier_gamma_calc(A_no, Gamma_N_m_TEV, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)
                
                # calculating terms or estimating first derivative
                
                F = np.array([abs(A[0]) - LESP, Gamma_tot_0])
                J = np.array([[(A_LEV_p[0] - A_LEV_m[0]) / (2*h), (A_TEV_p[0] - A_TEV_m[0]) / (4*2*h)],
                              [(Gamma_tot_p_LEV - Gamma_tot_m_LEV)/(2*h), (Gamma_tot_p_TEV - Gamma_tot_m_TEV)/(4*2*h)]])
                
                J_inv = np.linalg.inv(J)

                [Gamma_N[-1], Gamma_N[-2]] = np.array([x_i, y_i]) - J_inv@F 

        if t>0:
            print(A)
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

            # Finding induced velocity at each vortex
            for n in range(len(u_ind)):
                
                # Induced velocity on a vortex by another vortex
                for m in range(len(u_ind)):

                    u_ind_p, v_ind_p = pot.V_ind_ub(x_N[n], y_N[n], x_N[m], y_N[m],  Gamma_N[m], v_core)

                    u_ind[n] += u_ind_p
                    v_ind[n] += v_ind_p

                    if m == n:
                        u_ind[n] += 0.0
                        v_ind[n] += 0.0

                # Induced velocity on a vortex by the bounded vortex sheet            
                trans = lambda xi: np.arccos(1 - 2*xi/c)
                gamma = lambda xi: 2* U_ref * (A[0] * (1 + np.cos(trans(xi)))/np.sin(trans(xi)) + A[1] * np.sin(trans(xi)))# + A[2] * np.sin(2*trans(xi)) + A[3] * np.sin(3*trans(xi)) #+ A[4] * np.sin(4*trans(xi)) + A[5] * np.sin(5*trans(xi))

                u_ind_p, v_ind_p = pot.V_ind_b_fast_2(gamma, xi_N[n], eta_N[n], c, v_core)

                u_ind[n] += u_ind_p #+ U_ref
                v_ind[n] += v_ind_p #+ pot.hdot(t) 

            # Advecting the vortex
            x_N     = x_N + u_ind*t_step 
            y_N     = y_N + v_ind*t_step 

            # Shedding new TEV
            if t == t_step:
                x_N = np.append(x_N,(x_N[0] - (c-U_ref*t_step))*0.33 + c-U_ref*t_step)
                y_N = np.append(y_N,(y_N[0] - pot.h(t))*0.33 + pot.h(t))
            else:
                if LESP_flag:
                    x_N = np.append(x_N,(x_N[-2] - (c-U_ref*t))*0.33 + c-U_ref*t)
                    y_N = np.append(y_N,(y_N[-2] - pot.h(t))*0.33 + pot.h(t))
                else:
                    x_N = np.append(x_N,(x_N[-1] - (c-U_ref*t))*0.33 + c-U_ref*t)
                    y_N = np.append(y_N, pot.h(t))#(y_N[-1] - pot.h(t))*0.33 + pot.h(t))


            # Adjusting the body frame coordinates
            xi_N    = pot.xin2body(x_N, t, U_ref)
            eta_N   = pot.yin2body(y_N, t)

            # New TEV circulation
            Gamma_N = np.append(Gamma_N, 0.0)

            N += 1

            LESP_flag = 0

        # Calculate lift coefficient
        cl = np.append(cl, np.pi * (2 * A[0]+ A[1]))


        # # Movie
        # fig, ax = plt.subplots()
        # fig.dpi = 300
        # fig.set_size_inches(19.20, 10.80)
        # ax.plot(x_N, y_N, color='red', marker='o', markersize=10)
        # # ax.plot(xi_N, eta_N, 'bo')
        # # ax.plot([0, c], [0, 0], 'k')
        # ax.plot([0.0-U_ref *(t), c-U_ref*(t)], [pot.h(t), pot.h(t)], 'k')
        # ax.axis("equal")
        # ax.set_xlim(-30,5)
        # ax.set_ylim(-10,10)
        # buf = io.BytesIO()
        # plt.savefig(buf, format='png')
        # buf.seek(0)
        # frames.append(buf)
        # plt.close(fig)
    end = time.time()
    print(end - start)


    # Movie
    # images = []
    # for item in frames:
    #     image = plt.imread(item)
    #     images.append(image)

    # fig, ax = plt.subplots()

    # def update(frame):
    #     ax.clear()  # Clear the previous frame
    #     ax.imshow(images[frame], animated=True)
    #     ax.axis('off')
    #     return []
    
    # ani = animation.FuncAnimation(fig, update, len(images), blit=True)
    # writer = FFMpegWriter(fps=int(1/t_step))
    # ani.save('wtf.mp4', writer=writer,dpi=300)


    # print(t_step*U_ref/c)

    # plt.plot(t_d,cl)
    # plt.show()

    fig, ax = plt.subplots()
    fig.dpi = 300
    fig.set_size_inches(19.20, 10.80)
    ax.plot(x_N, y_N, 'ro')
    ax.plot([0.0-U_ref *(t_end), c-U_ref*(t_end)], [pot.h(t_end), pot.h(t_end)], 'k')
    ax.axis("equal")
    plt.show()

if __name__ == "__main__":
    main()




