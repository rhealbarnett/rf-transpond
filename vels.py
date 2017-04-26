#---------------------------------------------#
# drift velocities 			      #
# D Van Eester 2015 			      #
# "a crude model to study radio frequency     #
# induced density modification close to       #
# launchers"				      #
# figure 13, eqn A3                           #
# rlb 130417 				      #
#---------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt 
import pdb
from mpl_toolkits.mplot3d import Axes3D

T_ev = 5.
T = T_ev*11605.
B_mag = 3.4
me = 9.11e-31
qe = 1.6022e-19
kb = 1.38e-23
vt = np.sqrt(kb*T_ev / me)
om_c = qe*B_mag / me

#-- electric field
Ex = 0.
Ey = 0.
Ez = 0.
E = np.array([Ex, Ey, Ez])

#-- 'rotation' from e parallel to (ex, ey, ez); radians
alpha = 0.5
beta = 0.5

ex = np.sin(beta)*np.cos(alpha)
ey = np.sin(beta)*np.sin(alpha)
ez = np.cos(beta)

#ex = 0.0
#ey = 0.0
#ez = 1.0
e = np.array([ex, ey, ez])

b = e*B_mag

#-- initial velocity
v0 = vt/3.

vx = v0
vy = v0
vz = v0 
v = np.array([vx, vy, vz])

#-- initial position
xx = 0.0
xy = 0.0
xz = 0.0
x = np.array([xx, xy, xz])

t = 0.0
tmax = 1.0e-10 
dt = 1.0e-14

#---------------------------------------#
# RK4 -- unstable for large number of   #
# time steps.		                #
#---------------------------------------#

n = 2
m = 3

X = np.zeros((m, n))
X[:,0] = x
X[:,1] = v
print X

ax = 0.0
ay = 0.0
az = 0.0

a = np.array([0.0, 0.0, 0.0])

nmax = int((tmax - t)/dt)

def RHS(t, X):
	
	output = np.zeros((m,n))
	output[:,0] = X[:,1]
	output[:,1] = a + np.cross(X[:,1], b)

	return output

rkX_arr = np.zeros((nmax, m, n))

for ii in range(nmax):

	rkX_arr[ii,:,:] = X

	k1 = RHS(t, X)
	k2 = RHS(t + (dt/2.0), X + (dt/2.0)*k1)
	k3 = RHS(t + (dt/2.0), X + (dt/2.0)*k2)
	k4 = RHS(t + dt, X + dt*k3)

#	pdb.set_trace()

	X += (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)
	t += dt	

#	ax = 10.e6*(8.*X[0,0]**3.)
#	ay = 10.e6*(20.*X[1,0]**3.)
#	az = 10.e6*(40.*X[2,0]**3.)
#
#	a = np.array([ax, ay, az])

t_arr = np.linspace(0, tmax-dt, nmax)

#--------------------------------------#
# boris push method                    #
#--------------------------------------#

aE = (qe / me)*(dt / 2.0)*a
OM = (qe / me)*(dt / 2.0)*b
OM_mag_sq = OM[0]**2. + OM[1]**2. + OM[2]**2.
s = 2.*OM / (1.0 + OM_mag_sq)

boX_arr = np.zeros((nmax, m, n))

for ii in range(nmax):

	boX_arr[ii,:,0] = x
	boX_arr[ii,:,1] = v
	
	ax = 10.e6*(8.*x[0]**3.)
	ay = 10.e6*(20.*x[1]**3.)
	az = 10.e6*(40.*x[2]**3.)
	
	a = np.array([ax, ay, az])
	
	aE = (qe / me)*(dt / 2.0)*a
	
	vm = v + aE 
#	vp = vm + s*np.cross(v, OM) 
	vpr = vm + np.cross(vm, OM) 
	vpl = vm + np.cross(vpr, s) 
	v = vpl + aE 

	x = x - v*dt


#fig = plt.figure(1)
#ax = fig.gca(projection='3d')
#ax.plot(boX_arr[:,0,0], boX_arr[:,1,0], boX_arr[:,2,0])
#ax.set_title('Position (B0 = Bz, E0 = Ey, v0 = (vx,vy,vz))')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#plt.show(1)

plt.figure(2)
plt.plot(boX_arr[:,0,0], boX_arr[:,0,1], label='vx')
plt.plot(boX_arr[:,0,0], boX_arr[:,1,1], label='vy')
plt.plot(boX_arr[:,0,0], boX_arr[:,2,1], label='vz')
plt.xlabel('Position (x)')
plt.ylabel('Velocity (ms$^{-1}$)')
plt.legend(loc='upper right', prop={'size':10})
plt.show(2)




































