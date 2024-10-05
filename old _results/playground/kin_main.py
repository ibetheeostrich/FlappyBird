import math as mt
import numpy as np
import sympy as sym
import os

import subprocess as sp

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import helper as hp
import wing_kinematics as wk

import time



# Initialise kinematics solver ############################################

U = 10.0

# initialise root kinematics

rpm = 2

no_conf = 41

wing_lim = 15
x = np.linspace(0,4*np.pi,no_conf)
wing_angle_a = np.deg2rad(15)*np.cos(np.linspace(0,4*np.pi,no_conf))

wing_angle_w = -rpm*np.deg2rad(15)*np.sin(np.linspace(0,4*np.pi,no_conf))

# Variable and Array initialisation #######################################################

from geom import *

ainc_t = np.zeros(no_conf)
ainc_w = np.zeros(no_conf)

# initialise joint arrays

sh_a  = np.zeros((no_conf,3))  # shoulder joint

sb1_a = np.zeros((no_conf,3))  # shoulder blade joint
sb2_a = np.zeros((no_conf,3))  # shoulder blade joint

e1_a  = np.zeros((no_conf,3))  # elbow joint 1
e2_a  = np.zeros((no_conf,3))  # elbow joint 2
e3_a  = np.zeros((no_conf,3))  # elbow joint 3

w0_a  = np.zeros((no_conf,3))  # wrist root
w1_a  = np.zeros((no_conf,3))  # wrist joint 1
w2_a  = np.zeros((no_conf,3))  # wrist joint 2
w3_a  = np.zeros((no_conf,3))  # wrist joint 3
wt_a  = np.zeros((no_conf,3))  # wrist tip



# Establish range of motion numerically ############################################

motion = np.linspace(0,humerus_1_cam + scapula,100)
rom = []

for i in range(100):

    sba = np.array([0, motion[i], 0])

    try:
        config = wk.config_solver(sba)
        rom.append(motion[i])
    except:
        if len(rom) > 0:
            break
        
# Generating wing configs within range of motion as well as aoa at each chord section
# this section accounts for the camming used to keep the wing open on the down stroke
        
sba_a = np.zeros((no_conf,3)) # shoulder blade anchor

quart_index = int(np.floor(no_conf*0.25))
middle_index = int(np.floor(no_conf*0.5))
iiiquart_index = int(np.floor(no_conf*0.75))

sba_a[0:quart_index,1] = (max(rom) - min(rom))*0.5*np.cos(0) - (max(rom) + min(rom))*0.5
sba_a[middle_index:iiiquart_index,1] = (max(rom) - min(rom))*0.5*np.cos(0) - (max(rom) + min(rom))*0.5
sba_a[quart_index:middle_index,1] =  (max(rom) - min(rom))*0.5*np.cos(np.linspace(0,2*np.pi,quart_index )) - (max(rom) + min(rom))*0.5
sba_a[iiiquart_index:no_conf,1] =  (max(rom) - min(rom))*0.5*np.cos(np.linspace(0,2*np.pi,quart_index +1)) - (max(rom) + min(rom))*0.5


i = 0
for sba in sba_a:

    config = wk.config_solver(sba)

    # populating arrays
    sb1_a[i,:] = config[0]
    sb2_a[i,:] = config[1]

    e1_a[i,:]  = config[2]
    e2_a[i,:]  = config[3]
    e3_a[i,:]  = config[4]

    w0_a[i,:]  = config[5]
    w1_a[i,:]  = config[6]
    w2_a[i,:]  = config[7]
    w3_a[i,:]  = config[8]
    wt_a[i,:]  = config[9]

    # calculating wing tip velocity 
    V_t     = wing_angle_w[i]*wt_a[i,0]
    V_w     = wing_angle_w[i]*(np.sqrt(w0_a[i,0]**2 + w0_a[i,1]**2))

    # calculating relative angle of attack
    ainc_t[i]   = -np.arctan2(V_t,U)
    ainc_w[i] = -np.arctan2(V_w,U)

    # calculate planform area of wing at full extension
    if i == 0:
        s_ref = 0.5*((w0_a[i,0]*(root_chord + phalanx_3)) + 
                     ((wt_a[i,0] - w0_a[i,0])*(tip_chord + phalanx_3)))

    i += 1

# calculate wing configurations while "flapping"
    
# initialising array for intermediate  storage
conf = np.zeros((3,12))

for i in range(0,no_conf):
    # rearrange data for efficiency
    conf[:,0]  = sba_a[i,:]
    conf[:,1]  = sh_a[i,:]
    conf[:,2]  = sb1_a[i,:]
    conf[:,3]  = sb2_a[i,:]
    conf[:,4]  = e1_a[i,:]
    conf[:,5]  = e2_a[i,:]
    conf[:,6]  = e3_a[i,:]
    conf[:,7]  = w0_a[i,:]
    conf[:,8]  = w1_a[i,:]
    conf[:,9]  = w2_a[i,:]
    conf[:,10] = w3_a[i,:]
    conf[:,11] = wt_a[i,:]

    rot = np.array([
        [np.cos(wing_angle_a[i]), 0, -np.sin(wing_angle_a[i])],
        [0, 1, 0],
        [np.sin(wing_angle_a[i]), 0, np.cos(wing_angle_a[i])]
    ])

    new_conf = rot @ conf

    sba_a[i,:] = new_conf[:,0]  
    sh_a[i,:]  = new_conf[:,1]  
    sb1_a[i,:] = new_conf[:,2]  
    sb2_a[i,:] = new_conf[:,3]  
    e1_a[i,:]  = new_conf[:,4]  
    e2_a[i,:]  = new_conf[:,5]  
    e3_a[i,:]  = new_conf[:,6]  
    w0_a[i,:]  = new_conf[:,7]  
    w1_a[i,:]  = new_conf[:,8]  
    w2_a[i,:]  = new_conf[:,9]  
    w3_a[i,:]  = new_conf[:,10] 
    wt_a[i,:]  = new_conf[:,11] 

CL = np.zeros(int(np.ceil(no_conf/2)))

# Generating avl geometry files and running
for i in range(0,int(np.ceil(no_conf/2))):
    f = open("wing_template.txt", "r")
    template = f.readlines()
    f.close()

    template[6] = template[6].replace("SREF",str(s_ref))

    template[36] = template[36].replace("W0Y",str(w0_a[i,0]))
    template[36] = template[36].replace("W0X",str(-w0_a[i,1]))
    template[36] = template[36].replace("W0Z",str(w0_a[i,2]))
    template[36] = template[36].replace("AINCW",str(np.rad2deg(ainc_w[i])))

    template[43] = template[43].replace("WTY",str(wt_a[i,0]))
    template[43] = template[43].replace("WTX",str(-wt_a[i,1]))
    template[43] = template[43].replace("WTZ",str(wt_a[i,2]))
    template[43] = template[43].replace("AINCT",str(np.rad2deg(ainc_t[i])))

    # template[56] = template[56].replace("WTY",str(wt_a[i,0]))
    # template[56] = template[56].replace("WTX",str(-wt_a[i,1]))

    try:
        f = open("test{}.avl".format(i), "x")
    except:
        f = open("test{}.avl".format(i), "w")

    for line in template:
        f.write(line)
    f.close()


    # Running AVL for current config
    # initialise avl
    avl = sp.Popen(['/home/vikkiboi/FlappyBird/avl'],
               stdin=sp.PIPE,stdout=None, 
               stderr=None, 
               universal_newlines=True)

    # run AVL
    avl.stdin.write('LOAD test{}'.format(i,i) +'\n')
    avl.stdin.write('OPER\n' )
    avl.stdin.write('X\n')
    avl.stdin.write('S\n')
    avl.stdin.write('test{}.run'.format(i) +'\n')
    avl.stdin.write('\n')

    # close AVL
    avl.communicate('Quit\n')

    # delete geom file
    os.remove('test{}.avl'.format(i))

    # extract CL
    f = open("test{}.run".format(i), "r")
    bones = f.readlines()
    CL[i] = float(bones[15].split()[2])
    f.close()

    # delete results file
    os.remove('test{}.run'.format(i))

print('\n')

t = np.linspace(0,1/rpm,int(np.ceil(no_conf/2)))

plt.plot(t,CL)
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("C_L")
plt.show()






















    # plt.plot([sh_a[i,0], e1_a[i,0]],[sh_a[i,1], e1_a[i,1]])
    # plt.plot([sba_a[i,0], sb1_a[i,0]], [sba_a[i,1], sb1_a[i,1]])

    # plt.plot([sb2_a[i,0], e3_a[i,0]], [sb2_a[i,1], e3_a[i,1]])
    # plt.plot([e2_a[i,0], w0_a[i,0]], [e2_a[i,1], w0_a[i,1]])
    # plt.plot([e3_a[i,0], w1_a[i,0]], [e3_a[i,1], w1_a[i,1]])
    # plt.plot([w0_a[i,0], wt_a[i,0]], [w0_a[i,1], wt_a[i,1]])

    # ax = plt.gca()
    # ax.set_aspect('equal', adjustable='box')
    # ax.set_xlim(-0.5,4)
    # ax.set_ylim(-4,1)

    # plt.show()

    # panel preview

# def update(i):

#     ax.clear()
#     ax.set_aspect('equal')

#     ax.set_xlim(-6,6)
#     ax.set_ylim(-4,1)
#     ax.set_zlim(-4,4)

#     p1 = sh_a[i,:]
#     p2 = sh_a[i,:] + np.array([0, -1.0,0])
#     p3 = w0_a[i,:] + np.array([0, -0.5,0])
#     p4 = w0_a[i,:]
#     p5 = wt_a[i,:]
#     p6 = wt_a[i,:] + np.array([0, -0.1,0])

#     ax.plot([p1[0],p2[0]], [p1[1],p2[1]], [p1[2],p2[2]],'r')
#     ax.plot([p2[0],p3[0]], [p2[1],p3[1]], [p2[2],p3[2]],'r')
#     ax.plot([p3[0],p4[0]], [p3[1],p4[1]], [p3[2],p4[2]],'r')
#     ax.plot([p4[0],p1[0]], [p4[1],p1[1]], [p4[2],p1[2]],'r')

#     ax.plot([p4[0],p5[0]], [p4[1],p5[1]], [p4[2],p5[2]],'r')
#     ax.plot([p5[0],p6[0]], [p5[1],p6[1]], [p5[2],p6[2]],'r')
#     ax.plot([p6[0],p3[0]], [p6[1],p3[1]], [p6[2],p3[2]],'r')

#     ax.plot([-p1[0],-p2[0]], [p1[1],p2[1]], [p1[2],p2[2]],'r')
#     ax.plot([-p2[0],-p3[0]], [p2[1],p3[1]], [p2[2],p3[2]],'r')
#     ax.plot([-p3[0],-p4[0]], [p3[1],p4[1]], [p3[2],p4[2]],'r')
#     ax.plot([-p4[0],-p1[0]], [p4[1],p1[1]], [p4[2],p1[2]],'r')

#     ax.plot([-p4[0],-p5[0]], [p4[1],p5[1]], [p4[2],p5[2]],'r')
#     ax.plot([-p5[0],-p6[0]], [p5[1],p6[1]], [p5[2],p6[2]],'r')
#     ax.plot([-p6[0],-p3[0]], [p6[1],p3[1]], [p6[2],p3[2]],'r')

#     plt.tight_layout()

# # Create an animation
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.view_init(elev=0, azim=90, roll=0)

# ani = FuncAnimation(fig, update, np.arange(no_conf), interval = 0)

# plt.show()