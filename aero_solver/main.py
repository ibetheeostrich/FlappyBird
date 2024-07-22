import potential_flow as pot
import numpy as np

# define solver parameters

t_step = 0.01

t_span = np.arange(0.0,1.0,t_step)

# define free stream

u_inf = 10

v_inf = 0

# define camber line

no_points = 20

x_camber_line = np.linspace(0.0,1.0,no_points,True)

y_camber_line = np.zeros(no_points)

theta_camber_line = np.zeros(no_points)

# define wing motion

wing_u = np.sin(t_span*np.pi)

wing_v = np.zeros(len(t_span))



# solver

# initialise vortex positions
vort_x = x_camber_line
vort_y = y_camber_line

for t in t_span:

    if t == 0.0:

        X_n, X_m = np.meshgrid(vort_x,vort_x)

        Y_n, Y_m = np.meshgrid(vort_y,vort_y)

        theta_n, theta_m = np.meshgrid(theta_camber_line, theta_camber_line)

        int1 = np.multiply(X_n - X_m, X_n - X_m) + np.multiply(Y_n - Y_m, Y_n - Y_m) + np.identity(len(vort_x))

        int2 = np.multiply(X_n, np.cos(theta_m)) - np.multiply(Y_n, np.sin(theta_m)) - np.diag(np.diag(np.multiply(X_n, np.cos(theta_m)) - np.multiply(Y_n, np.sin(theta_m))))

        A1 = np.divide(int2,int1)

    break


