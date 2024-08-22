import helper as hp
import numpy as np
from geom import *

def config_solver(sba):

    # calulate sb1 position
    result = hp.circle_intersect(0.0, 0.0, shoulder_cam, 
                                 0.0, sba[1], scapula)

    for roots in result:

        if roots[0] > 0:
            sb1[0] = roots[0]
            sb1[1] = roots[1]
            break
        else:
            sb1[0] = roots[0]
            sb1[1] = roots[1]
            

    # calculate sb2 position
    sb2 = (sb1 - sba)/np.linalg.norm(sb1-sba)*scapula_cam + sba

    # calulate e1 position
    e1 = humerus_1*np.array([sb1[0],sb1[1],0])/np.linalg.norm(np.array([sb1[0],sb1[1],0]))

    # calculate e2 position
    a,b = hp.circle_intersect(e1[0], e1[1], humerus_1_cam, 
                                 sb2[0], sb2[1], humerus_2_cam)

    if a[1] < b[1]:
        e2[0] = a[0]
        e2[1] = a[1]
    else:
        e2[0] = b[0]
        e2[1] = b[1]

    # calculate e3
    e3 = (e2 - sb2)/np.linalg.norm((e2 - sb2))*humerus_2 + sb2

    # calculate w0
    w0 = (e1-e2)/np.linalg.norm(e1-e2)*radius + e2

    # calculate w1
    a,b = hp.circle_intersect(e3[0], e3[1], ulna, 
                                 w0[0], w0[1], phalanx_1_cam)

    if a[0] > b[0]:
        w1[0] = a[0]
        w1[1] = a[1]
    else:
        w1[0] = b[0]
        w1[1] = b[1]

    # calculate wt
    wt_int = (w1-w0)/np.linalg.norm(w1-w0)*phalanx_1 

    wt_int_xy = np.sqrt(np.linalg.norm(wt_int)**2 - wing_droop**2)
    
    wt = (w1-w0)/np.linalg.norm(w1-w0)*wt_int_xy + w0
    wt[2] = - wing_droop

    config = [sb1, sb2, e1, e2, e3, w0, w1, w2, w3, wt]

    flatten_angle = -np.arctan2(-wing_droop,wt[0])

    rot = np.array([
        [np.cos(flatten_angle), 0, -np.sin(flatten_angle)],
        [0, 1, 0],
        [np.sin(flatten_angle), 0, np.cos(flatten_angle)]
    ])  

    i=0
    for array in config:
        config[i] = rot @ array
        i += 1

    return config