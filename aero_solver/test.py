import math
import scipy
import scipy.integrate

def f():

    return lambda x,y: x**2 + y**1

def g(x):
    return x

def fg():

    a = f()

    return lambda x: a(x,3)

a = fg()

print(a(1))