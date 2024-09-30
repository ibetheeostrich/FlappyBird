import numpy as np
import scipy.integrate as inte


# Question 3
n7=2
n8=2
n9=6
A = 3 * (n8 + n9)

y_prime = lambda x: A/200 * (
    3 * (0.5*(1-np.cos(x)))**2 -
    4 * 0.5*(1-np.cos(x)) + 
    1
)

test = lambda x: A/200 * (
    3 * (0.5*(1-np.cos(x)))**2 -
    4 * 0.5*(1-np.cos(x)) + 
    1
) * (np.cos(x) - 1)

A1_func = lambda x: A/200 * (
    3 * (0.5*(1-np.cos(x)))**2 -
    4 * 0.5*(1-np.cos(x)) + 
    1
) * np.cos(x)

A2_func = lambda x: A/200 * (
    3 * (0.5*(1-np.cos(x)))**2 -
    4 * 0.5*(1-np.cos(x)) + 
    1
) * np.cos(2*x)


A1   = 2 / np.pi * inte.quad(A1_func,0,np.pi)[0]
A2   = 2 / np.pi * inte.quad(A2_func,0,np.pi)[0]

A0_req = -0.5 * A1

alpha0 = A0_req + inte.quad(y_prime,0,np.pi)[0] / np.pi 

print('Zero Lift AOA (degs):' + str(np.rad2deg(alpha0)))
print('Zero Lift AOA (rads):' + str(alpha0))

A0 = np.deg2rad(n7) - 1 / np.pi * inte.quad(y_prime,0,np.pi)[0]

CL = np.pi*(2*A0 + A1)
CM = 0.25*np.pi*(A2 - A1)

print('Coefficient of Lift:' + str(CL))
print('Moment Coefficient at 1/4:' + str(CM))

print(A0)
print(A1)
print(A2)
print(-1 / np.pi * inte.quad(test,0,np.pi)[0])