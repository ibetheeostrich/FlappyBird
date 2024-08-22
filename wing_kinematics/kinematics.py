import numpy as np
import helper as hp


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

class wing_kinematics:

    def __init__(self, I_in, I_out, II_in, II_out, III_in, III_out, IV_in, IV_out, V, VI_in, VI_out):

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

    def kin_2d(self, t, root_kin):
        '''
        t = current time step
        root_kin = kinematics of the root shoulder joint
        '''
        sa = np.array([0, 0])

        s1 = np.array([root_kin(t), 0])

        # find s3 position
        x1, x2, y1, y2 = hp.circle_intersect(sa[0], sa[1], self.III_in, s1[0], s1[1], self.I)

        if y1 > 0:
            s3 = np.array([x1, y1])
        elif y2 > 0:
            s3 = np.array([x2, y2])

        # find s2 position
        s2 = (s3 - s1) / np.linalg.norm(s3 - s1) * self.I_in

        # find e2
        e2 = (s3 - sa) / np.linalg.norm(s3 - sa) * self.III

        # find e1
        x1, x2, y1, y2 = hp.circle_intersect(e2[0], e2[1], self.IV_in, s2[0], s2[1], self.II_in)

        if x1 < x2:
            e1 = np.array([x1, y1])
        elif x1 > x2:
            e1 = np.array([x2, y2])

        # find e3
        e3 = (e1 - s2) / np.linalg.norm(e1 - s2) * self.II + s2

        # find w1
        w1 = (e2 - e1) / np.linalg.norm(e2 - e1) * self.IV + e1

        # find w2
        x1, x2, y1, y2 = hp.circle_intersect(e3[0], e3[1], self.V, w1[0], w1[1], self.VI_in)

        if y1 < y2:
            w2 = np.array([x2, y2])
        elif y1 > y2:
            w2 = np.array([x1, y1])

        # finding f3
        f3 = (w2 - w1) / np.linalg.norm(w2 - w1) * self.VI + w1

        return sa, s1, s2, s3, e1, e2, e3, w1, w2, f3