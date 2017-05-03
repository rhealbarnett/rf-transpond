#-----------------------------------#
# python rk4 package example        #
# from scipy.integrate		    #
# rlb 030517 			    #
#-----------------------------------#

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

T_ev = 5.0
T = T_ev*11605.
B_mag = 1.0
me = 9.11e-31
qe = 1.6022e-19
kb = 1.38e-23
vt = np.sqrt(kb*T_ev / me)
om_c = qe*B_mag / me

ex = 0.0
ey = 0.0
ez = 1.0

e = np.array([ex, ey, ez])

b = e*B_mag

xx = 0.0
xy = 0.0
xz = 0.0

x = np.array([xx, xy, xz])

v0 = vt/3.0

vx = v0
vy = v0
vz = v0

v = np.array([vx, vy, vz])

n = 2
m = 3

X = np.zeros((m, n))
X[:,0] = x
X[:,1] = v



def fun(x, z, a):
	"""
	Right hand side of the differential equations
	  dv/dt = a + (v cross b)
	  dx/dt = v
	"""
	x, v = z
	f = [v, a + np.cross(v,b)]
	return f

# create an ode instance to solve the system of differential equations defined
# by `fun', and set the solver method; dopri5 for rk4, dop853 for 8th order rk
solver = ode(fun)
solver.set_integrator('dopri5')

# give the value of omega to the solver. This is passed to `fun' when the 
# solver calls it
ax = 0.0
ay = 0.0
az = 0.0

a = np.array([ax, ay, az])
solver.set_f_params(a)

# set initial value z(0) = z0
t0 = 0.0
z0 = X
solver.set_initial_value(z0,t0)

# create the array `t' of time values at which to compute the solution, and create
# an array to hold the solution. Put the initial value in the solution array.
t1 = 1.0e-14
tmax = 1.0e-13
N = int((tmax - t0)/t1)
t = np.linspace(t0, tmax, N)
sol = np.zeros((N,m,n))
sol[0,:,:] = z0

# repeatedly call the `integrate' method to advance the solution to time t[k],
# and save the solution in sol[k]
k = 1
while solver.successful() and solver.t < tmax:
	solver.integrate(t[k])
	sol[k,:,:] = solver.y
	k += 1

# plot the solution
plt.plot(sol[:,0,0], sol[:,0,1], label = 'x')
#plt.plot(t, sol[:,1], label = 'y')
#plt.xlabel('t')
plt.legend()
plt.show()
