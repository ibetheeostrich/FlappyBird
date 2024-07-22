import math

def cart_vortex(x,y):

    u = 0.5 / math.pi() / (x**2 + y**2) * y 
    v = 0.5 / math.pi() / (x**2 + y**2) * x

    return u, v



    
