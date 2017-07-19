#-----------------------------------#
# accelerations toy problem         #
# D Van Eester 2015 		    #
# figure 12			    #
# rlb 290517  			    #
#-----------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import pdb

T_ev = 5.
B_mag = 3.4
u = 1.660539040e-27
m = 3.0160293*u
e = 1.6022e-19
vt = np.sqrt((T_ev*e) / m)
q = 2.*e
om_c = (q*B_mag) / m
period = (2.*np.pi) / om_c

#-- rotation matrix
alpha = 0.5
beta = 0.5

R = np.array([[np.cos(beta)*np.cos(alpha), np.cos(beta)*np.sin(alpha), -1.*np.sin(beta)],
		[-1.*np.sin(alpha), np.cos(alpha), 0.0],
		[np.sin(beta)*np.cos(alpha), np.sin(beta)*np.sin(alpha), np.cos(beta)]])

e_perp1 = R[0,:]
e_perp2 = R[1,:]
e_para = R[2,:]
b = B_mag*e_para

#-- initial velocity
v0 = vt/3.

#-- same v0 in x, y and z
vx = v0
vy = v0
vz = v0

v = np.array([vx, vy, vz])

#-- initial positions
x = 0.0
y = 0.0
z = 0.0

pos = np.array([x, y, z])

#-- acceleration
ax = 1.e6*(8.*pos[0]**3.)
ay = 1.e6*(20.*pos[1]**3.)
az = 1.e6*(40.*pos[2]**3.)

a = np.array([ax, ay, az])

#-- time step/iterations
dt = 1.0e-8
num_cyc = 5
num_ppc = 10
tmax = num_cyc*num_ppc*dt
#nmax = int(tmax/dt)
nmax = 18
num_dim = len(pos)
t = 0.0

#-- initialise arrays to store each value
a_arr = np.zeros((nmax, num_dim))
vd_arr = np.zeros((nmax, num_dim))
vcom_arr = np.zeros((nmax, num_dim))
pos_arr = np.zeros((nmax, num_dim))
X = np.zeros((num_dim, 2))
sol_arr = np.zeros((nmax, num_dim, 2))

X[:,0] = pos
X[:,1] = v

def RHS(t, pos):

	output = np.zeros((num_dim, 2))
	output[:,0] = X[:,1]
	output[:,1] = a 

	return output

for ii in range(nmax):

	vd = (1.0/om_c)*(np.cross(a, e_para))
	vcom = np.cross(vd, om_c*e_para)
	
	pos_arr[ii,:] = pos
	a_arr[ii,:] = a
	vd_arr[ii,:] = vd
	vcom_arr[ii,:] = vcom

#	k1 = RHS(t, X)
#	k2 = RHS(t + (dt/2.), X + (dt/2.)*k1)
#	k3 = RHS(t + (dt/2.), X + (dt/2.)*k2)
#	k4 = RHS(t + dt, X + dt*k3)
#
#	X += (dt/6.)*(k1 + 2.*k2 + 2.*k3 + k4)
#	t += dt
#
#	pos = X[:,0]
	ax = 1.e6*(8.*pos[0]**3.)
	ay = 1.e6*(20.*pos[1]**3.)
	az = 1.e6*(40.*pos[2]**3.)

	a = np.array([ax, ay, az])

	aE = dt*a

	v = v + aE
	pos = pos + v*dt

dv_dt = (vd_arr[1:nmax,:] - vd_arr[0:nmax-1,:])/dt

a_perp1 = np.sum(e_perp1*a_arr, axis=1)
a_perp2 = np.sum(e_perp2*a_arr, axis=1)
a_para = np.sum(e_para*a_arr, axis=1)
#vd_perp1 = (1.0/om_c)*(e_perp1*a_perp2)
#vd_perp2 = (1.0/om_c)*(e_perp2*a_perp1)


