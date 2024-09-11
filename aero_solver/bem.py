import math
import aero_solver.pot_func_simp as aero
import numpy as np
from copy import deepcopy
import scipy.integrate as inte


import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import io

# Constants and Globals

PI_inv = 1 / math.pi

def bem(tag,U_ref, alpha_eff, c, t_step, no_steps, kin):

    frames = []

    zeroth = np.array([])

    # Initialise problem
    t_d = np.arange(0,t_step*no_steps,t_step)
    cl = np.array([])

    pot = aero.aero_solver_osc_flat(kin, U_ref, t_step, alpha_eff)

    rho = 1.225

    lesp_flag = 0

    # Initialise vortex blob parameters
    no_gamma = 0
    Gamma_N = np.array([0.0])
    xi_N    = np.array([c[0]])
    eta_N   = np.array([0.0])
    x_N     = np.array([c[0]])
    y_N     = np.array([0.0])

    # Initialise Fourier coefficient matrix and calculation variabls
    A_no = 4
    # A = np.zeros(A_no)

    # Newton - Raphson Params
    dh = 0.01

    # Time loop

    for t in t_d:

        index = round(t/t_step)

        lesp = 0.2*c[index]

        # TEV Shedding
        if t > 0:

            if t == t_step:
                fourier_old = np.zeros(A_no)

            else:
                fourier_old = deepcopy(fourier)

            # Solving for TEV vorticity
            Gamma_err = 10000
  
            while abs(Gamma_err) > 0.00001:

                fourier, Gamma_sum, Gamma_tot = pot.fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, no_gamma, c[index], t)

                # print(Gamma_err, Gamma_b, sum(Gamma_N))
                Gamma_err = Gamma_tot

                if t > 0:

                    # Newton - Raphson iteration

                    # Guess
                    x_i = deepcopy(Gamma_N[-1])

                    # inputs at guess +h and -h to estimate first derivative
                    # calculating terms for estimating first derivative
                    Gamma_N_p = Gamma_N_m = Gamma_N

                    Gamma_N_p[-1] = x_i + dh
                    extra, Gamma_sum, Gamma_tot_p = pot.fourier_gamma_calc(A_no, Gamma_N_p, eta_N, xi_N, no_gamma, c[index], t)

                    Gamma_N_m[-1] = x_i - dh
                    extra, Gamma_sum, Gamma_tot_m = pot.fourier_gamma_calc(A_no, Gamma_N_m, eta_N, xi_N, no_gamma, c[index], t)

                    # Newton - Raphson iteration
                    Gamma_N[-1] = x_i - Gamma_tot / (0.5 * (Gamma_tot_p - Gamma_tot_m)/dh)

            # LEV Shedding
            if abs(fourier[0]) > lesp:

                if fourier[0] < 0:
                    lesp_c = - lesp

                else:
                    lesp_c = lesp

                stab = 0.1
                lesp_flag = 1
                Gamma_err = 100000

                Gamma_N = np.append(Gamma_N, -fourier[0]*10)

                # Gamma_N = np.append(Gamma_N, 10)
                # Gamma_N[-2] = 10

                xi_N = np.append(xi_N, 0)
                eta_N = np.append(eta_N, 0)

                x_N = np.append(x_N, pot.bodyin2x(xi_N[-1], t))
                y_N = np.append(y_N, pot.bodyin2y(eta_N[-1], t))

                no_gamma += 1

                iter_count = 0

                while abs(Gamma_err) > 0.00001 or abs(abs(fourier[0]) - lesp) > 0.00000001 :

                    fourier, Gamma_sum, Gamma_tot_0 = pot.fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, no_gamma, c[index], t)

                    # 2D Newton - Raphson iteration
                    # Guess
                    x_i = Gamma_N[-1] # LEV
                    y_i = Gamma_N[-2] # TEV

                    # inputs at guess +h and -h to estimate first derivative
                    Gamma_N_p_LEV = Gamma_N_m_LEV = Gamma_N_p_TEV = Gamma_N_m_TEV = Gamma_N

                    Gamma_N_p_LEV = np.array([y_i, x_i + stab * dh])
                    A_LEV_p, Gamma_sum, Gamma_tot_p_LEV = pot.fourier_gamma_calc(A_no, Gamma_N_p_LEV, eta_N[-2:], xi_N[-2:], 2, c[index], t)

                    Gamma_N_m_LEV = np.array([y_i, x_i - stab * dh])
                    A_LEV_m, Gamma_sum, Gamma_tot_m_LEV = pot.fourier_gamma_calc(A_no, Gamma_N_m_LEV, eta_N[-2:], xi_N[-2:], 2, c[index], t)

                    Gamma_N_p_TEV = np.array([y_i + stab * dh, x_i])
                    A_TEV_p, Gamma_sum, Gamma_tot_p_TEV = pot.fourier_gamma_calc(A_no, Gamma_N_p_TEV, eta_N[-2:], xi_N[-2:], 2, c[index], t)

                    Gamma_N_m_TEV = np.array([y_i - stab * dh, x_i])
                    A_TEV_m, Gamma_sum, Gamma_tot_m_TEV = pot.fourier_gamma_calc(A_no, Gamma_N_m_TEV, eta_N[-2:], xi_N[-2:], 2, c[index], t)

                    # calculating terms or estimating first derivative

                    F = np.array([fourier[0] - lesp_c, Gamma_tot_0])

                    J = np.array([[(A_LEV_p[0] - A_LEV_m[0]) / (2* stab *dh),         (A_TEV_p[0] - A_TEV_m[0]) / (2* stab *dh)],
                                  [(Gamma_tot_p_LEV - Gamma_tot_m_LEV)/(2* stab *dh), (Gamma_tot_p_TEV - Gamma_tot_m_TEV)/(2* stab *dh)]])

                    try:

                        # print(Gamma_N[-1])

                        J_inv = np.linalg.inv(J)

                        [Gamma_N[-1], Gamma_N[-2]] = np.array([x_i, y_i]) - J_inv@F 

                        Gamma_err = Gamma_tot_0

                        if iter_count > 100:
                            break
                        else:

                            iter_count += 1

                        # print(fourier)

                    except:

                        print(lesp, fourier)

                        print(J)

                        return cl, t_d[0:no_gamma-2], x_N, y_N, Gamma_N
                    
                            
        # Advecting and shedding vortices for next time step
        if t == 0:
            # xi_N    = np.array([c[index] + U_ref*t_step, c[index] + U_ref*t_step/3])
            xi_N    = np.array([c[index] + kin.pos_dot(t)*t_step, c[index] + kin.pos_dot(t)*t_step/3])
            eta_N   = np.array([kin.h(t_step), kin.h(t_step)/3])

            # x_N = np.array([c[index] , c[index] - U_ref*t_step/3])
            x_N = np.array([c[index] , c[index] - kin.pos_dot(t)*t_step/3])
            y_N = np.array([0, -kin.h(t_step)/3])

            Gamma_N = np.append(Gamma_N, 0.0)

            no_gamma += 2        

        if t > 0:

            # Calculating induced velocity on each vortex
            # u_ind, v_ind = pot.V_ind_tot_field(x_N,y_N,x_N,y_N, Gamma_N,fourier,no_gamma, U_ref,c[index],t)
            u_ind, v_ind = pot.V_ind_tot_field(x_N,y_N,x_N,y_N, Gamma_N,fourier,no_gamma, kin.pos_dot(t),c[index],t)

            # Advecting the vortex
            x_N     = x_N + u_ind*t_step 
            y_N     = y_N + v_ind*t_step 

            # # Shedding new TEV
            # if t == t_step:
            #     x_N = np.append(x_N,(x_N[0] - (c[index]-U_ref*t_step))*0.33 + c[index]-U_ref*t_step)
            #     y_N = np.append(y_N,(y_N[0] - kin.h(t))*0.33 + kin.h(t))
            # else:
            #     if lesp_flag:
            #         x_N = np.append(x_N,(x_N[-2] - (c[index]-U_ref*t))*0.33 + c[index]-U_ref*t)
            #         y_N = np.append(y_N,(y_N[-2] - kin.h(t))*0.33 + kin.h(t))
            #     else:
            #         x_N = np.append(x_N,(x_N[-1] - (c[index]-U_ref*t))*0.33 + c[index]-U_ref*t)
            #         y_N = np.append(y_N,(y_N[-1] - kin.h(t))*0.33 + kin.h(t))

            # Shedding new TEV
            if t == t_step:
                x_N = np.append(x_N,(x_N[0] - (c[index]-kin.pos(t)))*0.33 + c[index]-kin.pos(t))
                y_N = np.append(y_N,(y_N[0] - kin.h(t))*0.33 + kin.h(t))
            else:
                if lesp_flag:
                    x_N = np.append(x_N,(x_N[-2] - (c[index]-kin.pos(t)))*0.33 + c[index]-kin.pos(t))
                    y_N = np.append(y_N,(y_N[-2] - kin.h(t))*0.33 + kin.h(t))
                else:
                    x_N = np.append(x_N,(x_N[-1] - (c[index]-kin.pos(t)))*0.33 + c[index]-kin.pos(t))
                    y_N = np.append(y_N,(y_N[-1] - kin.h(t))*0.33 + kin.h(t))

            # Adjusting the body frame coordinates
            xi_N    = pot.xin2body(x_N, t)
            eta_N   = pot.yin2body(y_N, t)

            # New TEV circulation
            Gamma_N = np.append(Gamma_N, 0.0)

            no_gamma += 1

            lesp_flag = 0

        # Calculate lift coefficient

        # if t > 0:

            # disc_chord = np.linspace(0.0001,c[index], 500,endpoint = True)
            # disc_y = np.zeros(500)

            # lift_u, lift_v = pot.V_ind_ub_field(disc_chord, disc_y, xi_N, eta_N, Gamma_N, no_gamma)

            # disc_dphidx = lift_u*np.cos(alpha_eff) - lift_v*np.sin(alpha_eff)

            # trans = lambda xi: np.arccos(1 - 2*xi/c[index])
            # gamma = lambda xi: 2* U_ref * (fourier[0] * (1 + np.cos(trans(xi)))/np.sin(trans(xi)) + fourier[1] * np.sin(trans(xi)))# + fourier[2] * np.sin(2*trans(xi)) + fourier[3] * np.sin(3*trans(xi)) #+ fourier[4] * np.sin(4*trans(xi)) + fourier[5] * np.sin(5*trans(xi))

            # ub_terms = inte.trapezoid(disc_dphidx * gamma(disc_chord),disc_chord)

            # fourier_dot = (fourier - fourier_old)/t_step

            # f_n = rho * np.pi * c[index] * U_ref * (
            #     U_ref * np.cos(alpha_eff) * (fourier[0] + 0.5*fourier[1]) +
            #     c[index] * (0.75 * fourier_dot[0] + 0.25 * fourier_dot[1] + 0.125 * fourier_dot[2])
            # ) + rho * (
            #     ub_terms
            # )



            # cl = np.append(cl, f_n)

            index_close = np.where(xi_N < 1.5*c[index])

            cl_gamma = -np.pi * c[index] * U_ref * (2 * fourier[0]+ fourier[1]) #+ np.sum(Gamma_N[index_close])

            cl = np.append(cl, cl_gamma*1.225*U_ref*np.cos(alpha_eff))


        if t>0:
            x = np.linspace(0.0001, c, 513, endpoint=True)
            
            # cl = np.append(cl, np.pi * (2 * fourier[0]+ fourier[1]))
            # cl = np.append(cl, np.pi * (2 * fourier[0]+ fourier[1]))
            zeroth = np.append(zeroth, fourier[0])

            print(t, t_step/c[index]*U_ref,  np.pi * c[index] * U_ref * (fourier[0] + fourier[1] * 0.5))


            # Pressure Field
            # x = np.linspace(-5.2,0.5,200)
            # y = np.linspace(-1.5,1.5,100)

            # X,Y = np.meshgrid(x,y)

            # X_straight = np.reshape(X,-1)
            # Y_straight = np.reshape(Y,-1)

            # U,V = pot.V_ind_tot_field(X_straight, Y_straight, x_N, y_N, Gamma_N,fourier,no_gamma, U_ref,c[index],t)

            # U = np.reshape(U,newshape=(100,200))
            # V = np.reshape(V,newshape=(100,200))

            # cp = - (U**2 + V**2) / U_ref**2

            # fig, ax = plt.subplots()
            # fig.dpi = 300
            # fig.set_size_inches(19.20, 10.80)
            # contf = ax.contourf(X,Y,cp,levels=np.linspace(-1.0, 1.0, 100), extend='both')
            # fig.colorbar(contf,
            #            orientation='horizontal',
            #            shrink=0.5, pad = 0.1,
            #            ticks=[-1.0, -1.0, 0.0, 1.0])
            # # ax.plot(x_N, y_N, 'ro')
            # ax.plot([0.0-kin.pos(t), c[index]-kin.pos(t)], [kin.h(t), kin.h(t)], 'k')
            # ax.axis("equal")
            # ax.set_xlim(-5.2,0.5)
            # ax.set_ylim(-1.5,1.5)
            # plt.savefig(str(index) + 'pressure'+'.png',)
            # plt.close(fig)

            # Movie
            fig, ax = plt.subplots()
            fig.dpi = 300
            fig.set_size_inches(19.20, 10.80)
            ax.plot(x_N, y_N, 'ro')
            # ax.plot(xi_N, eta_N, 'bo')
            # ax.plot([0, c], [0, 0], 'k')
            ax.plot([0.0-kin.pos(t), c[index]-kin.pos(t)], [kin.h(t), kin.h(t)], 'k')
            ax.axis("equal")
            # ax.set_xlim(-20.2,0.5)
            # ax.set_ylim(-1.5,1.5)
            plt.savefig(str(index) + '.png',)
            plt.close(fig)


    return tag, cl, t_d[0:-1], x_N, y_N, Gamma_N, zeroth
    






