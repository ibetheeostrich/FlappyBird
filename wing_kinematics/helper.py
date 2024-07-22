import numpy as np
import math as mt

def vec_rot(v: np.ndarray,k: np.ndarray,theta: float):
    '''
    VEC_ROT - computes the rotation of vector v around unit vector k about the angle
    theta according to the right hand rule
    '''
    # calculate rotation
    v_rot = v * np.cos(theta) + np.cross(k, v) * np.sin(theta) + k * np.dot(k, v) * (1 - np.cos(theta))

    return v_rot

def circle_intersect(x0, y0, r0, x1, y1, r1):
    '''
    CIRCLE_INTERSECT: calculates the intersection of two circles

    source: https://stackoverflow.com/questions/55816902/finding-the-intersection-of-two-circles
    '''

    # circle 1: (x0, y0), radius r0
    # circle 2: (x1, y1), radius r1

    d=mt.sqrt((x1-x0)**2 + (y1-y0)**2)
    
    # non intersecting
    if d > r0 + r1 :
        return None
    # One circle within other
    if d < abs(r0-r1):
        return None
    # coincident circles
    if d == 0 and r0 == r1:
        return None
    else:
        a=(r0**2-r1**2+d**2)/(2*d)
        h=mt.sqrt(r0**2-a**2)
        x2=x0+a*(x1-x0)/d   
        y2=y0+a*(y1-y0)/d   
        x3=x2+h*(y1-y0)/d     
        y3=y2-h*(x1-x0)/d 

        x4=x2-h*(y1-y0)/d
        y4=y2+h*(x1-x0)/d
        
        return (x3, y3), (x4, y4)
    
    