import math
import aero_solver.pot_func as aero
import numpy as np
from copy import deepcopy
import scipy.integrate as inte

# Constants and Globals

PI_inv = 1 / math.pi

def bem(U_ref, alpha_eff, c, t_step, no_steps, kin):

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
    A_no = 3
    # A = np.zeros(A_no)

    # Newton - Raphson Params
    dh = 0.01

    # Time loop

    for t in t_d:

        index = int(t/t_step)

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

                stab = -1

                lesp_flag = 1
                Gamma_err = 100000

                Gamma_end = deepcopy(Gamma_N[-1])

                Gamma_N = np.append(Gamma_N, fourier[0]*0.5)
                # Gamma_N[-2] = fourier[0]

                # Gamma_N = np.append(Gamma_N, 10)
                # Gamma_N[-2] = 10

                xi_N = np.append(xi_N, 0)
                eta_N = np.append(eta_N, - y_N[-1]/abs(y_N[-1]) * 0.01 * c[index] )

                x_N = np.append(x_N, pot.bodyin2x(xi_N[-1], t))
                y_N = np.append(y_N, pot.bodyin2y(eta_N[-1], t))

                no_gamma += 1

                while abs(Gamma_err) > 0.00001 or abs(abs(fourier[0]) - lesp) > 0.001 :

      

                    fourier, Gamma_sum, Gamma_tot_0 = pot.fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, no_gamma, c[index], t)

                    # 2D Newton - Raphson iteration
                    # Guess
                    x_i = deepcopy(Gamma_N[-1])
                    y_i = deepcopy(Gamma_N[-2])

                    # inputs at guess +h and -h to estimate first derivative
                    Gamma_N_p_LEV = Gamma_N_m_LEV = Gamma_N_p_TEV = Gamma_N_m_TEV = Gamma_N

                    Gamma_N_p_LEV[-1] = x_i + stab * dh
                    Gamma_N_p_LEV[-2] = y_i
                    A_LEV_p, Gamma_sum, Gamma_tot_p_LEV = pot.fourier_gamma_calc(A_no, Gamma_N_p_LEV, eta_N, xi_N, no_gamma, c[index], t)

                    Gamma_N_m_LEV[-1] = x_i - stab * dh
                    Gamma_N_m_LEV[-2] = y_i
                    A_LEV_m, Gamma_sum, Gamma_tot_m_LEV = pot.fourier_gamma_calc(A_no, Gamma_N_m_LEV, eta_N, xi_N, no_gamma, c[index], t)

                    Gamma_N_p_TEV[-2] = y_i + stab * dh
                    Gamma_N_p_TEV[-1] = x_i
                    A_TEV_p, Gamma_sum, Gamma_tot_p_TEV = pot.fourier_gamma_calc(A_no, Gamma_N_p_TEV, eta_N, xi_N, no_gamma, c[index], t)

                    Gamma_N_m_TEV[-2] = y_i - stab * dh
                    Gamma_N_m_TEV[-1] = x_i
                    A_TEV_m, Gamma_sum, Gamma_tot_m_TEV = pot.fourier_gamma_calc(A_no, Gamma_N_m_TEV, eta_N, xi_N, no_gamma, c[index], t)

                    # calculating terms or estimating first derivative

                    # F = np.array([abs(fourier[0]) - lesp, Gamma_tot_0])

                    F = np.array([abs(fourier[0]) - lesp, Gamma_tot_0])

                    J = np.array([[(A_LEV_p[0] + A_LEV_m[0]) / (2* stab *dh),         (A_TEV_p[0] + A_TEV_m[0]) / (2* stab *dh)],
                                  [(Gamma_tot_p_LEV - Gamma_tot_m_LEV)/(2* stab *dh), (Gamma_tot_p_TEV - Gamma_tot_m_TEV)/(2* stab *dh)]])

                    try:

                        J_inv = np.linalg.inv(J)

                        [Gamma_N[-1], Gamma_N[-2]] = np.array([x_i, y_i]) - J_inv@F 

                        Gamma_err = Gamma_tot_0

                        print(fourier)

                    except:

                        return cl, t_d[0:no_gamma-2], x_N, y_N
                            
        # Advecting and shedding vortices for next time step
        if t == 0:
            xi_N    = np.array([c[index] + U_ref*t_step, c[index] + U_ref*t_step/3])
            eta_N   = np.array([kin.h(t_step), kin.h(t_step)/3])

            x_N = np.array([c[index] , c[index] - U_ref*t_step/3])
            y_N = np.array([0, -kin.h(t_step)/3])

            Gamma_N = np.append(Gamma_N, 0.0)

            no_gamma += 2

        if t > 0:

            # Calculating induced velocity on each vortex
            u_ind, v_ind = pot.V_ind_tot_field(x_N,y_N,x_N,y_N, Gamma_N,fourier,no_gamma, U_ref,c[index],t)

            # Advecting the vortex
            x_N     = x_N + u_ind*t_step 
            y_N     = y_N + v_ind*t_step 

            # Shedding new TEV
            if t == t_step:
                x_N = np.append(x_N,(x_N[0] - (c[index]-U_ref*t_step))*0.33 + c[index]-U_ref*t_step)
                y_N = np.append(y_N,(y_N[0] - kin.h(t))*0.33 + kin.h(t))
            else:
                if lesp_flag:
                    x_N = np.append(x_N,(x_N[-2] - (c[index]-U_ref*t))*0.33 + c[index]-U_ref*t)
                    y_N = np.append(y_N,(y_N[-2] - kin.h(t))*0.33 + kin.h(t))
                else:
                    x_N = np.append(x_N,(x_N[-1] - (c[index]-U_ref*t))*0.33 + c[index]-U_ref*t)
                    y_N = np.append(y_N,(y_N[-2] - kin.h(t))*0.33 + kin.h(t))#(y_N[-1] - pot.h(t))*0.33 + pot.h(t))


            # Adjusting the body frame coordinates
            xi_N    = pot.xin2body(x_N, t)
            eta_N   = pot.yin2body(y_N, t)

            # New TEV circulation
            Gamma_N = np.append(Gamma_N, 0.0)

            no_gamma += 1

            lesp_flag = 0

        # Calculate lift coefficient

        if t > 0:

            disc_chord = np.linspace(0.0001,c[index], 500,endpoint = True)
            disc_y = np.zeros(500)

            lift_u, lift_v = pot.V_ind_ub_field(disc_chord, disc_y, xi_N, eta_N, Gamma_N, no_gamma)

            disc_dphidx = lift_u*np.cos(alpha_eff) - lift_v*np.sin(alpha_eff)

            trans = lambda xi: np.arccos(1 - 2*xi/c[index])
            gamma = lambda xi: 2* U_ref * (fourier[0] * (1 + np.cos(trans(xi)))/np.sin(trans(xi)) + fourier[1] * np.sin(trans(xi)))# + fourier[2] * np.sin(2*trans(xi)) + fourier[3] * np.sin(3*trans(xi)) #+ fourier[4] * np.sin(4*trans(xi)) + fourier[5] * np.sin(5*trans(xi))

            ub_terms = inte.trapz(disc_dphidx * gamma(disc_chord),disc_chord)

            fourier_dot = (fourier - fourier_old)/t_step

            f_n = rho * np.pi * c[index] * U_ref * (
                U_ref * np.cos(alpha_eff) * (fourier[0] + 0.5*fourier[1]) +
                c[index] * (0.75 * fourier_dot[0] + 0.25 * fourier_dot[1] + 0.125 * fourier_dot[2])
            ) + rho * (
                ub_terms
            )



            # cl = np.append(cl, np.pi * (2 * fourier[0]+ fourier[1]))
            cl = np.append(cl, f_n)
            # print(fourier[0])
        else:
            cl = np.append(cl,0)

        

    return cl, t_d, x_N, y_N

    






