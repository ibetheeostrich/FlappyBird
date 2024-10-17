import numpy as np
import scipy.integrate as inte

import matplotlib.pyplot as plt
import matplotlib

from aero_solver.bem import bem_2 as bem
from wing_kinematics.kinematics import blade_element_kinematics as bek
  
from wing_kinematics.kinematics import wing_kinematics as wk

from multiprocessing import Pool

from geom import *

import time

import os

import csv

# Initialise problem based off wind tunnel data

wt_path = './windtunnel_results/DYNAMIC ALL (CLEAN1P)/'
wt_results = os.listdir(wt_path)

case = wt_results[3]

params = case.split('_')

rho = 1.225
U_ref       = float(params[1][:-2])
alpha       = float(params[0][:-3])
alpha_eff   = np.deg2rad(alpha)    
frequency   = float(params[2][:-6])

alpha_eff   = np.deg2rad(0.0) 
U_ref       = 10
frequency   = 2.0

scale = 0.5

no_bem = 12 



print(params)
print(U_ref,alpha,frequency)

# Time span
no_steps = 400

t_step = 1/frequency/no_steps

t_span = np.linspace(0.0, no_steps*t_step, no_steps, endpoint=False)

amp = 42.5

lesp = 10.25

# get data from csv

wt_lift = [] 
wt_drag = [] 
wt_time = []

with open(wt_path + case,'r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    for row in lines: 
        wt_lift.append(row[1]) 
        wt_drag.append(row[0])
        wt_time.append(row[2])

wt_lift = np.array(wt_lift[1:], dtype=np.float32)  
wt_drag = np.array(wt_drag[1:], dtype=np.float32)  
wt_time = np.array(wt_time[1:], dtype=np.float32)

def wing_plot():
    t = 0.0
    wing_kin = wk(I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out, F_I, F_II, A_I, A_II)

    root_pos = -0.07

    points = wing_kin.kin_2d_plot(root_pos)

    print(points)

    np.savetxt("points1.csv", points,  
              delimiter = ",")  

def main():
    
    start = time.time()

    # Initialises the geometry and allows the kinematics to be solved
    wing_kin = wk(I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out, F_I, F_II, A_I, A_II)

    # Defining root kinematics that will drive morphing
    root_kin = lambda x: 1*(- 0.0585 - 0.0125*np.sin(2 * np.pi * frequency * x)) 


    # Initialising array for storing BEM kinematics
    r_pos_temp = np.zeros((no_bem,no_steps+1))  
    chords_temp = np.zeros((no_bem,no_steps+1))
    le_pos_temp = np.zeros((no_bem,no_steps+1))
    area_temp = np.zeros(no_steps+1)

    for t in np.append(t_span,t_span[-1]+t_step):

        sa, w1, f1, f2, f3 = wing_kin.kin_2d_V2(t, root_kin)


        # calculating the chord, blade element position and leading position at each time step
        r_pos, chords, le_pos, area = wing_kin.variable_wing_params(sa, w1, f1, f2, f3, no_bem)

        r_pos_temp[:, round(t/t_step)] = r_pos
        chords_temp[:, round(t/t_step)] = chords
        area_temp[round(t/t_step)] = area
        le_pos_temp[:, round(t/t_step)] = le_pos

    # Initialise input list
    args = []

    # Assigning tag to order the blade elements
    tag = np.linspace(0,no_bem,no_bem,endpoint=False)

    # Calculating leading edge positions
    for i in range(no_bem-1):

        le_pos_temp[i, :] = le_pos_temp[i, :] - le_pos_temp[i, 0]

    # Calculating BEM kinematics and populating argument array
    for i in range(no_bem-1):

        kin = bek(np.deg2rad(amp) , frequency, r_pos_temp[i,:], chords_temp[i,:], le_pos_temp[i,:], U_ref, alpha_eff, t_step)

        args.append((round(tag[i]),U_ref, alpha_eff, chords_temp[i,:], t_step, no_steps, kin, lesp))

    # Multiprocessing
    def pool_handler():
        p = Pool(16)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        results = p.starmap(bem, args)

        return results

    # Evaluating blade elements using multi processing
    a = pool_handler()

    # Results accumulation matrices 
    cl_mat = np.zeros((no_bem,no_steps))
    cd_mat = np.zeros((no_bem,no_steps))

    r_mat = r_pos_temp[:,:-1]

    # Processing results for each blade element
    for results in a:

        # Getting BEM tag
        i = results[0]

        # Formatting Data
        cl_mat[i,:] = cl = results[1]
        cd_mat[i,:] = cd = results[2]

        td = results[3]

        x = results[4]
        y = results[5]

        gamma = results[6]
        pos_arr = np.where(gamma > 0)[0]
        neg_arr = np.where(gamma < 0)[0]

        # Creating Figures
        fig, ax = plt.subplots()
        fig.dpi = 300
        ax.plot(td, cl)
        ax.plot(td, cd)
        ax.set_xlabel('Time  (s)')
        ax.set_ylabel('Lift Force (n)')
        plt.savefig('l v t' + str(i) + '.png')
        plt.close(fig)

    # Integrating BEM
    l_int = inte.trapezoid(cl_mat*np.cos(np.deg2rad(amp)*np.cos(2*np.pi*frequency*t_span)),r_mat,axis=0)
    d_int = inte.trapezoid(cd_mat,r_mat,axis=0)

    print(inte.trapezoid(l_int[2:],td[2:])) 
    print(inte.trapezoid(d_int[2:],td[2:])) 
     

    fig, ax = plt.subplots()
    fig.dpi = 300
    ax.plot(td[:-2], l_int[2:],'r')
    # ax.plot(td[:-2], -d_int[2:],'g')

    ax.plot(wt_time,wt_lift,'b')
    # ax.plot(wt_time,wt_drag,'y')

    ax.set_xlabel('Time  (s)')
    ax.set_ylabel('Lift Force (n)')
    plt.savefig('lint v t' + str(i) + '.png')
    plt.clf()

    fig, ax = plt.subplots()
    fig.dpi = 300
    ax.plot(td[:-2], d_int[2:],'r')
    # ax.plot(td[:-2], -d_int[2:],'g')

    # ax.plot(wt_time,wt_lift,'b')
    # ax.plot(wt_time,wt_drag,'y')

    ax.set_xlabel('Time  (s)')
    ax.set_ylabel('Drag Force (n)')
    plt.savefig('dint v t' + str(i) + '.png')
    plt.clf()

    print(time.time()-start)

    np.savetxt(f"{t_step:.7f}_{alpha:.1f}deg_{U_ref:.1f}ms_{frequency:.1f}Hz_{lesp:.2f}LESP" + '.csv', np.transpose(np.vstack((-d_int,l_int,td))),  
              delimiter = ",")

if __name__ == "__main__":
    main()
    # wing_plot()