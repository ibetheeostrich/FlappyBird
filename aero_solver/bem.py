import math
import aero_solver.pot_func_simp as aero
import numpy as np
from copy import deepcopy
import scipy.integrate as inte
import aero_solver.aero_objects as ao

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import io

# Constants and Globals

PI_inv = 1 / math.pi

def bem_old(tag,U_ref, alpha_eff, c, t_step, no_steps, kin):

    frames = []

    zeroth = np.array([])

    # Initialise problem
    t_d = np.arange(0,t_step*no_steps,t_step)
    cl = np.array([])

    pot = aero.aero_solver_osc_flat(kin, U_ref, t_step, alpha_eff,np.max(c))

    rho = 1.225

    lesp_flag = 0

    lev_shed_flag = 0

    # Initialise vortex blob parameters
    no_gamma = 0
    Gamma_N = np.array([0.0])
    xi_N    = np.array([c[0]])
    eta_N   = np.array([0.0])
    x_N     = np.array([c[0]])
    y_N     = np.array([0.0])

    # Initialise Fourier coefficient matrix and calculation variabls
    A_no = 35
    # A = np.zeros(A_no)

    # Newton - Raphson Params
    dh = 0.001

    # Time loop

    for t in t_d:

        index = round(t/t_step)

        lesp = 0.15

        # TEV Shedding
        if t > 0:

            if t == t_step:
                fourier_old = np.zeros(A_no)

            else:
                fourier_old = deepcopy(fourier)

            # Solving for TEV vorticity
            Gamma_err = 10000
  
            while abs(Gamma_err) > 0.00001:

                fourier, Gamma_sum, Gamma_tot = pot.fourier_gamma_calc_2(A_no, Gamma_N, eta_N, xi_N, no_gamma, c[index], t)

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
                    # extra, Gamma_sum, Gamma_tot_p = pot.fourier_gamma_calc_2(A_no, Gamma_N_p, eta_N, xi_N, no_gamma, c[index], t)
                    extra, Gamma_sum, Gamma_tot_p = pot.fourier_gamma_calc_2(A_no, [x_i + dh], [eta_N[-1]], [xi_N[-1]], 1, c[index], t)

                    Gamma_N_m[-1] = x_i - dh
                    # extra, Gamma_sum, Gamma_tot_m = pot.fourier_gamma_calc_2(A_no, Gamma_N_m, eta_N, xi_N, no_gamma, c[index], t)
                    extra, Gamma_sum, Gamma_tot_m = pot.fourier_gamma_calc_2(A_no, [x_i - dh], [eta_N[-1]], [xi_N[-1]], 1, c[index], t)

                    # Newton - Raphson iteration
                    Gamma_N[-1] = x_i - Gamma_tot / (0.5 * (Gamma_tot_p - Gamma_tot_m)/dh)

            if abs(fourier[0]) < lesp:
               lev_shed_flag = 0

            # LEV Shedding
            if abs(fourier[0]) > lesp:

                if fourier[0] < 0:
                    lesp_c = -lesp
                    Gamma_N = np.append(Gamma_N, 10)

                else:
                    lesp_c = lesp
                    Gamma_N = np.append(Gamma_N,-10)

                stab = 0.1
                lesp_flag = 1
                Gamma_err = 100000

                if lev_shed_flag:
                    xi_N = np.append(xi_N, xi_N[-2]*0.33 )
                    eta_N = np.append(eta_N, eta_N[-2]*0.33 )

                elif not lev_shed_flag:

                    if lesp < 0:

                        xi_N = np.append(xi_N, c[index]*0.005 )
                        eta_N = np.append(eta_N, c[index]*0.005 )

                    elif lesp > 0:

                        xi_N = np.append(xi_N, c[index]*0.005 )
                        eta_N = np.append(eta_N, -c[index]*0.005 )

                lev_shed_flag = 1

                x_N = np.append(x_N, pot.bodyin2x(xi_N[-1], t))
                y_N = np.append(y_N, pot.bodyin2y(eta_N[-1], t))

                no_gamma += 1

                iter_count = 0

                # while abs(Gamma_err) > 0.00001 or abs(abs(fourier[0]) - lesp) > 0.00000001 :
                while True:

                    fourier, Gamma_sum, Gamma_tot_0 = pot.fourier_gamma_calc_2(A_no, Gamma_N, eta_N, xi_N, no_gamma, c[index], t)

                    # 2D Newton - Raphson iteration
                    # Guess
                    x_i = deepcopy(Gamma_N[-1]) # LEV
                    y_i = deepcopy(Gamma_N[-2]) # TEV

                    Gamma_N_p_LEV = np.array([y_i, x_i + stab * dh])
                    A_LEV_p, Gamma_sum, Gamma_tot_p_LEV = pot.fourier_gamma_calc_2(A_no, Gamma_N_p_LEV, eta_N[-2:], xi_N[-2:], 2, c[index], t)

                    Gamma_N_m_LEV = np.array([y_i, x_i - stab * dh])
                    A_LEV_m, Gamma_sum, Gamma_tot_m_LEV = pot.fourier_gamma_calc_2(A_no, Gamma_N_m_LEV, eta_N[-2:], xi_N[-2:], 2, c[index], t)

                    Gamma_N_p_TEV = np.array([y_i + stab * dh, x_i])
                    A_TEV_p, Gamma_sum, Gamma_tot_p_TEV = pot.fourier_gamma_calc_2(A_no, Gamma_N_p_TEV, eta_N[-2:], xi_N[-2:], 2, c[index], t)

                    Gamma_N_m_TEV = np.array([y_i - stab * dh, x_i])
                    A_TEV_m, Gamma_sum, Gamma_tot_m_TEV = pot.fourier_gamma_calc_2(A_no, Gamma_N_m_TEV, eta_N[-2:], xi_N[-2:], 2, c[index], t)

                    # calculating terms or estimating first derivative

                    target = np.array([fourier[0] - lesp_c, Gamma_tot_0])

                    jacob = np.array([[(A_LEV_p[0] - A_LEV_m[0]) / (2* stab *dh),         (A_TEV_p[0] - A_TEV_m[0]) / (2* stab *dh)],
                                  [(Gamma_tot_p_LEV - Gamma_tot_m_LEV)/(2* stab *dh), (Gamma_tot_p_TEV - Gamma_tot_m_TEV)/(2* stab *dh)]])

                    try:

                        # print(Gamma_N[-1])

                        J_inv = np.linalg.inv(jacob)

                        [Gamma_N[-1], Gamma_N[-2]] = np.array([x_i, y_i]) - J_inv@target 

                        Gamma_err = Gamma_tot_0

                        if iter_count > 1000:
                            break
                        else:

                            iter_count += 1

                        if abs(Gamma_err) < 0.0000001 and abs(abs(fourier[0]) - lesp) < 0.00000001 and iter_count > 5:
                            break 

                        # print(fourier)

                    except:

                        print(lesp, fourier)

                        print(jacob)

                        return cl, t_d[0:no_gamma-2], x_N, y_N, Gamma_N
                    
                # print(Gamma_N[-1], fourier[0], Gamma_tot_0)
                    

        # Calculate lift coefficient

        if t > 0:

            fourier_dot = (fourier - fourier_old)/t_step




            theta = np.linspace(0,np.pi,70,endpoint=True)

            xi = 0.5*c[index]*(1-np.cos(theta))

            eta = np.zeros(70)

            u_ind, junk = pot.V_ind_ub_field(xi, eta, xi_N, eta_N, Gamma_N, 1)

            fourier_inf = fourier[0]*(1+np.cos(theta))

            for i in range(1,len(fourier)):

                fourier_inf += fourier[i]*np.sin(i*theta)*np.sin(theta)

            fourier_inf = 2*U_ref*np.cos(alpha_eff) * fourier_inf

            vort_cont = inte.trapezoid(0.5*c[index]*fourier_inf*(u_ind + kin.pos_dot(t)),theta) / (0.5 * c[index] * (U_ref*np.cos(alpha_eff))**2)            

            c_n = vort_cont + 2*np.pi*(
                (fourier[0] + 0.5*fourier[1]) + 
                c[index]/(U_ref) * (0.75 * fourier_dot[0] + 0.25 * fourier_dot[1] + 0.125 * fourier_dot[2])
            ) 


            cl = np.append(cl, c_n)

            index_close = np.where(xi_N < 1.5*c[index])



            if tag == 4:
                # print(fourier[0],fourier[1],fourier[2])
                # print(np.pi*(2*fourier[0] + fourier[1]), fourier)
                # print(F,iter_count)

                print(c_n,fourier[0],fourier[1],kin.h_dot(t))


        if t>0:
            x = np.linspace(0.0001, c, 513, endpoint=True)
            
            # cl = np.append(cl, np.pi * (2 * fourier[0]+ fourier[1]))
            # cl = np.append(cl, np.pi * (2 * fourier[0]+ fourier[1]))
            zeroth = np.append(zeroth, fourier[0])

            # print(t, t_step/c[index]*U_ref,  np.pi * c[index] * U_ref * (fourier[0] + fourier[1] * 0.5))


            # # Pressure Field
            if tag == 4:
                x = np.linspace(-0.1,0.5,100)
                y = np.linspace(-0.25,0.25,100)   

                X,Y = np.meshgrid(x,y)  

                X_straight = np.reshape(X,-1)
                Y_straight = np.reshape(Y,-1)   

                U,V = pot.V_ind_tot_field(X_straight, Y_straight, xi_N, eta_N, Gamma_N,fourier,no_gamma, U_ref,c[index],t) 


                U = np.reshape(U,newshape=(100,100)) + kin.pos_dot(t)
                V = np.reshape(V,newshape=(100,100)) - kin.h_dot(t)   

                cp =  1 - (U**2 + V**2) / ((U_ref)**2) 

                # cp =  U**2


                fig, ax = plt.subplots()
                fig.dpi = 300
                fig.set_size_inches(10.80, 10.80)
                contf = ax.contourf(X,Y,cp,levels=np.linspace(-2.0, 1.0, 100), extend='both')
                fig.colorbar(contf,
                           orientation='horizontal',
                           shrink=0.5, pad = 0.1,
                           ticks=[0.0, 1.0])
                ax.plot(xi_N, eta_N, 'ro', ms=1)
                ax.plot([0.0, c[index]], [0, 0], 'k')
                ax.axis("equal")
                ax.set_xlim(-0.1,0.5)
                ax.set_ylim(-0.1,0.1)
                plt.savefig('pressure'+str(index)+'.png',)
                plt.close(fig)

            # Pressure Field @ vort
            # if tag == 4:

            #     U,V = pot.V_ind_tot_field(xi_N, eta_N, xi_N, eta_N, Gamma_N,fourier,no_gamma, U_ref,c[index],t) 


            #     U += kin.pos_dot(t)
            #     V -= kin.h_dot(t)   

            #     cp =  1 - (U**2 + V**2) / ((U_ref)**2) 

            #     # cp =  U**2


            #     fig, ax = plt.subplots()
            #     fig.dpi = 300
            #     fig.set_size_inches(10.80, 10.80)
            #     contf = ax.scatter(xi_N,eta_N,c=cp)#,levels=np.linspace(-2.0, 1.0, 100), extend='both')
            #     fig.colorbar(contf,orientation='horizontal')
            #     ax.plot(xi_N, eta_N, 'ro', ms=1)
            #     ax.plot([0.0, c[index]], [0, 0], 'k')
            #     ax.axis("equal")
            #     ax.set_xlim(-0.1,0.5)
            #     ax.set_ylim(-0.1,0.1)
            #     plt.savefig('pressure @ vort'+str(index)+'.png',)
            #     plt.close(fig)

            # # Movie
            # if tag ==4:
            #     fig, ax = plt.subplots()
            #     fig.dpi = 300
            #     fig.set_size_inches(19.20, 10.80)
            #     contf = ax.scatter(x_N,y_N,c=Gamma_N)#,levels=np.linspace(-2.0, 1.0, 100), extend='both')
            #     fig.colorbar(contf,orientation='horizontal')
            #     # ax.plot(xi_N, eta_N, 'bo')
            #     # ax.plot([0, c], [0, 0], 'k')
            #     ax.plot([0.0-kin.pos(t), c[index]-kin.pos(t)], [kin.h(t), kin.h(t)], 'k')
            #     ax.axis("equal")
            #     # ax.set_xlim(-20.2,0.5)
            #     # ax.set_ylim(-1.5,1.5)
            #     plt.savefig('vorticity'+str(index)+'.png',)
            #     plt.close(fig)

            # Vorticity distribution
            # if tag ==4:
            #     fig, ax = plt.subplots()
            #     fig.dpi = 300
            #     fig.set_size_inches(19.20, 10.80)
            #     ax.plot(disc_chord, gamma(disc_chord))
            #     ax.set_xlim(0,c[index])
            #     ax.set_ylim(-20,20)
            #     plt.savefig('vorticity_dist'+str(index)+'.png',)
            #     plt.close(fig)                            
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
            u_ind, v_ind = pot.V_ind_tot_field(xi_N,eta_N,xi_N,eta_N, Gamma_N,fourier,no_gamma, kin.pos_dot(t),c[index],t)

            xi_N      = xi_N  + u_ind*t_step + kin.pos_dot(t)*t_step
            eta_N     = eta_N + v_ind*t_step - kin.h_dot(t)*t_step

            # Shedding new TEV
            if t == t_step:
                xi_N = np.append(xi_N,(xi_N[0] - c[index])*0.33 + c[index])
                eta_N = np.append(eta_N,(eta_N[0])*0.33)
            else:
                if lesp_flag:
                    xi_N = np.append(xi_N,(xi_N[-2] - c[index])*0.33 + c[index])
                    eta_N = np.append(eta_N,(eta_N[-2])*0.33)
                else:
                    xi_N = np.append(xi_N,(xi_N[-1] - c[index])*0.33 + c[index])
                    eta_N = np.append(eta_N,(eta_N[-1])*0.33)


            # # Adjusting the body frame coordinates
            # xi_N    = pot.xin2body(x_N, t)
            # eta_N   = pot.yin2body(y_N, t)

            # Adjusting the body frame coordinates
            x_N    = pot.bodyin2x(xi_N, t)
            y_N    = pot.bodyin2y(eta_N, t)

            # New TEV circulation
            Gamma_N = np.append(Gamma_N, 0.0)

            no_gamma += 1

            lesp_flag = 0


    return tag, cl, t_d[0:-1], x_N, y_N, Gamma_N, zeroth
    

def bem(tag,U_ref, alpha_eff, chords, t_step, no_steps, kin, lesp_crit):

    cl = np.zeros(no_steps)
    cd = np.zeros(no_steps)

    td = np.linspace(0,no_steps*t_step,no_steps,endpoint=False)

    x_dot = lambda t: U_ref - (kin.le[round(t/t_step) + 1] - kin.le[round(t/t_step)]) / t_step
    h_dot = lambda t: 2 * np.pi * kin.f * kin.r[round((t) / t_step)] * kin.aa * np.sin(2 * np.pi * kin.f * (t)) #if t > 0.02 else 0.0
    alpha_dot = lambda t: 0.0 

    u = lambda t: U_ref*t - kin.le[round(t/t_step)] 
    h = lambda t: kin.r[round((t) / t_step)] * kin.aa - kin.r[round((t) / t_step)] * kin.aa * np.cos(2 * np.pi * kin.f * (t)) #if t > 0.02 else 0.0
    alpha = lambda t: alpha_eff

    be  = ao.camber_line(chords, 35, x_dot,h_dot,alpha_dot,u,h,alpha,t_step)

    field   = ao.vorticity_field(chords[0])

    lev_flag = 0

    for t in td:

        if t > -1:

            lev_flag_prev = lev_flag

            be.fourier_old = deepcopy(be.fourier) 

            be.update_pos(t)

            field.shed_tev(be)

            be.kelvinkutta(field,0.001,t)

            if abs(be.fourier[0]) > lesp_crit:          

                field.shed_lev(be)

                be.kelvinlesp(field, 0.001, lesp_crit, t)

                gb = np.pi * be.c(t) * be.x_dot(t) * (be.fourier[0] + be.fourier[1] * 0.5) + np.sum(np.concatenate((field.tev, field.lev, field.ext)))           
                print(f'{t:.2f} {abs(be.fourier[0]) - lesp_crit:.3f} {gb}') 

                lev_flag = 1

            if lev_flag_prev == 1 and lev_flag == 0:

                field.ext   = np.append(field.ext   ,deepcopy(field.lev)    )
                field.ext_x = np.append(field.ext_x ,deepcopy(field.lev_x)  )
                field.ext_y = np.append(field.ext_y ,deepcopy(field.lev_y)  )

                field.lev   = np.array([])
                field.lev_x = np.array([])
                field.lev_y = np.array([])
    #####################################################################################    

            # if round(t/t_step) % 5 == 0 and tag == 9:
            #     fig, ax = plt.subplots()
            #     fig.dpi = 300
            #     ax.scatter(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
            #             np.concatenate((field.tev_y, field.lev_y, field.ext_y))
            #             , c='b', s=1.7)
            #     ax.plot(be.x,
            #             be.y,
            #             'k')
            #     ax.axis("equal")
            #     ax.set_xlim(be.x[0] - 0.1,be.x[-1] + 0.1)
            #     ax.set_ylim(be.y[0] - 0.1,be.y[-1] + 0.1)
            #     plt.savefig(str(round(t/t_step)) + '.png')
            #     plt.clf()   


            # print(x_dot(t))

    #####################################################################################    

            field.advect(be,t_step,t)

            cl[round(t/t_step)], cd[round(t/t_step)] = be.calc_cl(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                                             np.concatenate((field.tev_y, field.lev_y, field.ext_y)),
                                             np.concatenate((field.tev, field.lev, field.ext)),
                                             t, t_step)


    return tag, cl, cd, td, np.concatenate((field.tev_x, field.lev_x, field.ext_x)), np.concatenate((field.tev_y, field.lev_y, field.ext_y)), np.concatenate((field.tev, field.lev, field.ext)), 1


def bem_span(tag, chords, t_step, no_steps, freq, amp):

    cl = np.zeros(no_steps)
    cd = np.zeros(no_steps)

    td = np.linspace(0,no_steps*t_step,no_steps,endpoint=False)

    x_dot = lambda t: 0.0000000000000000001
    h_dot = lambda t: 0.0
    alpha_dot = lambda t: 2 * np.pi * freq * amp * np.sin(2 * np.pi * freq * t)

    u = lambda t: 0.0
    h = lambda t: 0.0
    alpha = lambda t: amp - amp * np.cos(2 * np.pi * freq * t)

    lesp_crit = 100.15

    be  = ao.camber_line(chords, 35, x_dot,h_dot,alpha_dot,u,h,alpha,t_step)

    field   = ao.vorticity_field(chords[0])

    for t in td:

        if t > 0:

            be.fourier_old = deepcopy(be.fourier)   

            be.update_pos(t)

            field.shed_tev(be)

            be.kelvinkutta(field,0.001,t)

            # if abs(be.fourier[0]) > lesp_crit:          

            #     field.shed_lev(be)

            #     be.kelvinlesp(field, 0.001, lesp_crit, t)

    #####################################################################################    

            if True:
                fig, ax = plt.subplots()
                fig.dpi = 300
                ax.plot(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                        np.concatenate((field.tev_y, field.lev_y, field.ext_y))
                        ,'ro')
                ax.plot(be.x,
                        be.y,
                        'k')
                ax.axis("equal")
                plt.savefig(str(round(t/t_step)) + '.png')
                plt.clf()   

                            
            # print(t) 

    #####################################################################################    

            field.advect(be,t_step,t)

            cl[round(t/t_step)], cd[round(t/t_step)] = be.calc_cl(np.concatenate((field.tev_x, field.lev_x, field.ext_x)),
                                             np.concatenate((field.tev_y, field.lev_y, field.ext_y)),
                                             np.concatenate((field.tev, field.lev, field.ext)),
                                             t, t_step)


    return tag, cl, cd, td
