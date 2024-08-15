import math
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

                A[i] *= 2.0 / np.pi / U_ref

    Gamma_b = np.pi * c * U_ref * (A[0] + A[1] * 0.5)

    return A, sum(Gamma_N), Gamma_b + sum(Gamma_N)


A_no = 2
Gamma_N = eta_N = xi_N = np.array([])
c = 1

N = 0

v_core = 0.1

U_ref = 2
alpha_eff = 0.0


g_trans = lambda theta: 0.5 * c * (1 - np.cos(theta))

t = 0.04

A, B, C= fourier_gamma_calc(A_no, Gamma_N, eta_N, xi_N, U_ref, alpha_eff, v_core, g_trans, c, N, t)

print(C - B)