import math
import numpy as np
import scipy.integrate as inte
from multiprocessing import Pool
import time
import matplotlib.pyplot as plt
import pot_func as pot

## AI SLOP

# Constants and Globals
PI_INV = 1 / math.pi

def int_dphideta(params):
    v_core = 5 * 1.3 * 0.05 * 2
    alpha_eff = 0
    c = 2.0

    Gamma_n, eta_n, xi_n = params
    dphideta = pot.dphideta(xi_n, eta_n, Gamma_n, v_core, alpha_eff)
    
    g_trans = lambda theta: 0.5 * c * (1 - np.cos(theta))
    integrand_n = lambda theta: dphideta(g_trans(theta))
    
    A_int, _ = inte.quad(integrand_n, 0.0, np.pi)
    return A_int

def main():
    U_REF = 2
    U_REF_INV = 1 / U_REF   
    alpha_eff = np.deg2rad(0)   
    c = 2.0
    t_step = 0.05
    t_end = 1000 * t_step
    t_d = np.arange(0, t_end, t_step)
    cl = []

    N = 0
    Gamma_N = [0.0]
    xi_N = [c]
    eta_N = [0.0]
    x_N = [c]
    y_N = [0.0]
    v_core = 5 * 1.3 * t_step * U_REF

    g_trans = lambda theta: 0.5 * c * (1 - np.cos(theta))
    eta_xi = lambda xi: 0.0

    A = np.zeros(6)
    int_bounds = np.linspace(0, np.pi, num=500)

    for t in t_d:
        start = time.time()
        Gamma_err = 10000
  
        while abs(Gamma_err) > 0.001:
            W_0 = pot.W_0(U_REF, alpha_eff, t)

            for i in range(len(A)):
                A[i] = 0.0

                if i == 0:
                    if N == 0:
                        A[i] = -PI_INV * U_REF_INV * pot.W_0_fast_1(U_REF, alpha_eff, t) * np.pi
                    else:
                        A[i] = pot.W_0_fast_1(U_REF, alpha_eff, t) * np.pi
                        if N < 170:
                            for n in range(N):                            
                                Gamma_n = Gamma_N[n]
                                eta_n = eta_N[n]
                                xi_n = xi_N[n]
                                dphideta = pot.dphideta(xi_n, eta_n, Gamma_n, v_core, alpha_eff)
                                integrand_n = lambda theta: dphideta(g_trans(theta))
                                A_int, _ = inte.quad(integrand_n, 0.0, np.pi)
                                A[i] += A_int
                        else:
                            inter = np.array([Gamma_N, eta_N, xi_N]).T
                            with Pool() as pool:
                                tempcoeff = pool.map(int_dphideta, inter)
                            A[i] += sum(tempcoeff)
                    A[i] *= -1.0 / np.pi / U_REF
                else:
                    integrand_n = lambda theta: W_0(theta) * math.cos(i * theta)
                    A[i], _ = inte.quad(integrand_n, 0.0, np.pi)
                    A[i] *= 2.0 / np.pi / U_REF
                    if N > 0:
                        for n in range(N):
                            Gamma_n = Gamma_N[n]
                            eta_n = eta_N[n]
                            xi_n = xi_N[n]
                            dphideta = pot.dphideta(xi_n, eta_n, Gamma_n, v_core, alpha_eff)
                            integrand_n = lambda theta: dphideta(g_trans(theta)) * math.cos(i * theta)
                            A_int, _ = inte.quad(integrand_n, 0.0, np.pi)
                            A[i] += A_int
                        A[i] *= 2.0 / np.pi / U_REF
            
            Gamma_b = np.pi * c * U_REF * (A[0] + A[1] * 0.5)
            Gamma_err = Gamma_b + sum(Gamma_N)
            Gamma_N[-1] -= 0.5 * Gamma_err

        if t == 0:
            xi_N = [c + U_REF * t_step]
            eta_N = [pot.h(t_step)]
        if t > 0:
            u_ind = np.zeros(N+1)
            v_ind = np.zeros(N+1)
            for n in range(N+1):
                for m in range(N+1):
                    u_ind_p, v_ind_p = pot.V_ind_ub(x_N[n], y_N[n], x_N[m], y_N[m], Gamma_N[m], v_core)
                    u_ind[n] += u_ind_p
                    v_ind[n] += v_ind_p
            x_N = x_N + u_ind * t_step
            y_N = y_N + v_ind * t_step

            if t == t_step:
                x_N = np.append(x_N, (x_N[0] - (c - U_REF * t_step)) * 0.33 + c - U_REF * t_step)
                y_N = np.append(y_N, (y_N[0] - pot.h(t_step)) * 0.33 + pot.h(t_step))
            else:
                x_N = np.append(x_N, (x_N[-1] - (c - U_REF * t)) * 0.33 + c - U_REF * t)
                y_N = np.append(y_N, (y_N[-1] - pot.h(t)) * 0.33 + pot.h(t))

            xi_N = pot.xin2body(x_N, t, U_REF)
            eta_N = pot.yin2body(y_N, t)

            Gamma_N = np.append(Gamma_N, 0.0)
            N += 1

        end = time.time()
        cl.append(np.pi * (2 * A[0] + A[1]))
        V = pot.hdot(t)
        print(cl[-1], np.rad2deg(np.arctan2(V, U_REF)), A[0], N, end - start)

    plt.plot(x_N, y_N, 'ro')
    plt.plot([0.0 - U_REF * t_end, c - U_REF * t_end], [pot.h(t_end), pot.h(t_end)], 'k')
    plt.plot(xi_N, eta_N, 'bo')
    plt.axis("equal")
    plt.show()

    plt.plot(t_d, cl)
    plt.show()

    print(t_step * U_REF / c)

if __name__ == "__main__":
    main()
