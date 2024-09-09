import numpy as np
import math as mt
import shapely

# Define some parameters
'''
Nomenclature

s - shoulder joints
e - elbow joints
w - wrist joints
f - fingers

roman numerals: lengths

in arrays 
position [0] = x
position [1] = y
'''

sa = np.array([0, 0])

s1 = np.array([0, 0])
s2 = np.array([0, 0])
s3 = np.array([0, 0])

class blade_element_kinematics:
    def __init__(self, amplitude, freq, span_pos, chord, le_pos, U_ref, t_step):

        self.U_ref = U_ref

        self.aa = amplitude     # radians
        self.f = freq           # hz?
        self.r = span_pos       # m

        self.chord = chord      # m
        self.le = le_pos        # m 

        self.wavelength = 1.0 / freq

        self.t_step = t_step

    def h(self, t):

        index = round(t / self.t_step)

        a = self.r[index] * self.aa

        # if t < 0.25:
        #     return 0
    
        # if t >= 0.0 and t < 0.5 * self.wavelength:

        #     return 0.5*(a - a * np.cos(2 * np.pi * self.f * (t)))
        
        # elif t >= 0.5 * self.wavelength:
        #     return a * np.cos(2 * np.pi * self.f * (t - 0.5/self.f)) 

        # if t<0.25:
        #     return 0
        # else:
        #     return -7.0*(t - 0.25)

        return a * np.cos(2 * np.pi * self.f * (t - 0.25/self.f)) 

    def h_dot(self, t):

        index = round(t / self.t_step)

        a = self.r[index] * self.aa

        # if t < 0.25:
        #     return 0
        
        # if t >= 0.0 and t < 0.5 * self.wavelength:
        #     return np.pi * self.f * a * np.sin(2 * np.pi * self.f * (t))
        
        # elif t >= 0.5 * self.wavelength:
        #     return - 2 * np.pi * self.f * a * np.sin(2 * np.pi * self.f * (t - 0.5/self.f)) 
        # if t<0.25:
        #     return 0
        # else:
        #     return -7.0

        return - 2 * np.pi * self.f * a * np.sin(2 * np.pi * self.f * (t - 0.25/self.f)) 
        
#################################################################################
#               QUARANTINE                                                      #
#################################################################################

    def pos(self, t):

        index = round(t / self.t_step)

        return self.U_ref*t - self.le[index]
    
    def pos_dot(self, t):

        index = round(t / self.t_step)

        return self.U_ref - (self.le[index + 1] - self.le[index]) / self.t_step

    # def c(self, t):

    #     index = int(t / self.t_step)

    #     return self.chord[index]

    # def c_dot(self, t):

    #     index = int(t / self.t_step)

    #     return self.chord[index+1] - self.chord[index]
        
#################################################################################



class wing_kinematics:

    def __init__(self, I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out,F_I, F_II, A_I, A_II):

        self.I_in   = I_in
        self.II_in  = II_in
        self.III_in = III_in
        self.IV_in  = VI_in
        self.VI_in  = VI_in

        self.I_out   = I_out
        self.II_out  = II_out
        self.III_out = III_out
        self.IV_out  = VI_out
        self.VI_out  = VI_out

        self.I      = I_in + I_out
        self.II     = II_in + II_out
        self.III    = III_in + III_out
        self.IV     = IV_in + IV_out
        self.V      = V
        self.VI     = VI_in + VI_out

        self.F_I = F_I
        self.F_II = F_II

        self.A_I = A_I
        self.A_II = A_II

    def kin_2d_V2(self, t, root_kin):
        '''
        t = current time step
        root_kin = kinematics of the root shoulder joint
        '''

        # Shoulder positions
        sa = np.array([0, 0])

        s1 = np.array([root_kin(t), 0])

        s3 = np.array([
            - (self.I**2 - self.III_in**2 - s1[0]**2) * 0.5 / s1[0],
            np.sqrt(self.III_in**2 - ((self.I**2 - self.III_in**2 - s1[0]**2) * 0.5 / s1[0])**2)
        ])

        s2 = (s3 -s1) * self.I_in / self.I + s1


        # Elbow positions
        e2 = s3 * self.III / self.III_in

        temp_D = np.linalg.norm(e2-s2)

        temp_X = - (self.IV_in**2 - self.II_in**2 - temp_D**2) / (2*temp_D)

        temp_Y = np.sqrt(self.II_in**2 - temp_X**2)

        e1 = np.array([
            (temp_X * (e2[0] - s2[0]) + temp_Y * (s2[1] - e2[1])) / temp_D + s2[0],
            (temp_X * (e2[1] - s2[1]) + temp_Y * (e2[0] - s2[0])) / temp_D + s2[1]
        ])

        e3 = (e1 - s2) * self.II/self.II_in + s2

        # Look at the flick of the wrist
        w1 = (e2-e1) * self.IV/self.IV_in + e1

        temp_D = np.linalg.norm(w1-e3)

        temp_X = - (self.VI_in**2 - self.V**2 - temp_D**2) / (2*temp_D)

        temp_Y = np.sqrt(self.V**2 - temp_X**2)

        w2 = np.array([
            (temp_X * (w1[0] - e3[0]) + temp_Y * (e3[1] - w1[1])) / temp_D + e3[0],
            (temp_X * (w1[1] - e3[1]) + temp_Y * (w1[0] - e3[0])) / temp_D + e3[1]
        ])        

        # Fingertips
        f3 = (w2 - w1) * self.VI / self.VI_in + w1

        int1 = (w2 - e3) / np.linalg.norm(w2 - e3) * (self.V-self.A_I) + e3
        f1 = (int1 - w1) / np.linalg.norm(int1 - w1) * self.F_I + w1

        int1 = (w2 - e3) / np.linalg.norm(w2 - e3) * (self.V-self.A_II) + e3
        f2 = (int1 - w1) / np.linalg.norm(int1 - w1) * self.F_II + w1

        return sa, w1, f1, f2, f3


    def kin_2d(self, t, root_kin):
        '''
        t = current time step
        root_kin = kinematics of the root shoulder joint
        '''
        sa = np.array([0, 0])

        s1 = np.array([root_kin(t), 0])

        # find s3 position
        x1, y1, x2, y2 = circle_intersect(sa[0], sa[1], self.III_in, s1[0], s1[1], self.I)

        if y1 > 0:
            s3 = np.array([x1, y1])
        elif y2 > 0:
            s3 = np.array([x2, y2])

        # find s2 position
        s2 = (s3 - s1) / np.linalg.norm(s3 - s1) * self.I_in + s1

        # find e2
        e2 = (s3 - sa) / np.linalg.norm(s3 - sa) * self.III

        # find e1
        x1, y1, x2, y2 = circle_intersect(e2[0], e2[1], self.IV_in, s2[0], s2[1], self.II_in)

        if x1 < x2:
            e1 = np.array([x1, y1])
        elif x1 > x2:
            e1 = np.array([x2, y2])

        # find e3
        e3 = (e1 - s2) / np.linalg.norm(e1 - s2) * self.II + s2

        # find w1
        w1 = (e2 - e1) / np.linalg.norm(e2 - e1) * self.IV + e1

        # find w2
        x1, y1, x2, y2 = circle_intersect(e3[0], e3[1], self.V, w1[0], w1[1], self.VI_in)

        if y1 < y2:
            w2 = np.array([x2, y2])
        elif y1 > y2:
            w2 = np.array([x1, y1])

        # finding f3
        f3 = (w2 - w1) / np.linalg.norm(w2 - w1) * self.VI + w1

        # finding f1
        int1 = (w2 - e3) / np.linalg.norm(w2 - e3) * (self.V-self.A_I) + e3
        f1 = (int1 - w1) / np.linalg.norm(int1 - w1) * self.F_I + w1

        # finding f2
        int1 = (w2 - e3) / np.linalg.norm(w2 - e3) * (self.V-self.A_II) + e3
        f2 = (int1 - w1) / np.linalg.norm(int1 - w1) * self.F_II + w1

        return sa, w1, f1, f2, f3
    
    def variable_wing_params(self, sa, w1, f1, f2, f3, no_bem):
        '''
        calculates the chord as the wing morphs as well as the change in radial distance
        no_bem must be an even number
        '''

        

        wing_vert = [tuple(sa), tuple(w1), tuple(f3), tuple(f2), tuple(f1), (-0.08, 0.0),tuple(sa)]

        wing_poly = shapely.Polygon(wing_vert)

        r_pos1 = np.linspace(0.0, w1[1], int(no_bem*0.5)+1)

        r_pos2 = np.linspace(w1[1], f3[1], int(no_bem*0.5) + 1, endpoint=True)

        r_pos = np.hstack([r_pos1[1:], r_pos2[1:]])

        chords = np.zeros(len(r_pos))

        le_pos = np.zeros(len(r_pos))

        i = 0
        for r in r_pos:

            if r == r_pos[-1]:
                
                chords[i] = 0

                le_pos[i] = 0

            else:                

                line = [(20, r), (-20, r)]
                shapely_line = shapely.geometry.LineString(line)

                intersection = list(wing_poly.intersection(shapely_line).coords)

                chords[i] = abs(intersection[0][0] - intersection[1][0])

                le_pos[i] = max(intersection[0][0], intersection[1][0])

            i += 1

        return r_pos, chords, le_pos, wing_poly.area

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
        
        return x3, y3, x4, y4
    
    