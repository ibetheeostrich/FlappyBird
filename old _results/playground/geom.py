import numpy as np

# initialise control linkage lengths

scapula     = 0.4
humerus_1   = 1.0
humerus_2   = 0.8
radius      = 1.5
ulna        = 1.48
phalanx_1   = 2
phalanx_2   = 0.5
phalanx_3   = 0.5

# initialise cam joint positions

shoulder_cam    = 0.25
scapula_cam     = 0.2
humerus_1_cam   = 0.2
humerus_2_cam   = 0.6
phalanx_1_cam   = 0.1

# initialise other geometric constraints
wing_droop  = 0.2
root_chord  = 1.0
tip_chord   = 0.1

###############################################################################

# vector init
sh  = [0,0,0]   # shoulder joint

sb1 =np.array([0.0,0.0,0.0])   # shoulder blade joint
sb2 =np.array([0.0,0.0,0.0])   # shoulder blade joint

e1 = np.array([0.0,0.0,0.0])   # elbow joint 1
e2 = np.array([0.0,0.0,0.0])   # elbow joint 2
e3 = np.array([0.0,0.0,0.0])   # elbow joint 3

w0 = np.array([0.0,0.0,0.0])   # wrist root
w1 = np.array([0.0,0.0,0.0])   # wrist joint 1
w2 = np.array([0.0,0.0,0.0])   # wrist joint 2
w3 = np.array([0.0,0.0,0.0])   # wrist joint 3
wt = np.array([0.0,0.0,0.0])   # wrist tip
