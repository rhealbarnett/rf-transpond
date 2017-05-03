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
period = 1.0 / (om_c / (2.*np.pi))

#-- electric field
Ex = 0.
Ey = 0.
Ez = 0.
E = np.array([Ex, Ey, Ez])

#-- rotation matrix; radians
alpha = 0.5
beta = 0.5

R = np.array([[np.cos(beta)*np.cos(alpha), np.cos(beta)*np.sin(alpha), -1.*np.sin(beta)],
		[-1.*np.sin(alpha), np.cos(alpha), 0.0],
		[np.sin(beta)*np.cos(alpha), np.sin(beta)*np.sin(alpha), np.cos(beta)]])

e_para = R[2,:]

b = e_para*B_mag

#-- initial velocity
v0 = vt/3.

vx = v0 
vy = v0
vz = v0
v = np.array([vx, vy, vz])

rl = me*v[1]/(qe*b[2])

#-- initial position
xx = 0.0
xy = 0.0
xz = 0.0
x = np.array([xx, xy, xz])

t = 0.0
tmax = 1.0e-10 
dt = 1.0e-14

#---------------------------------------#
# RK4    			        #
#---------------------------------------#

n = 2
h = 3

X = np.zeros((h, n))
X[:,0] = x
X[:,1] = v

#X = np.zeros(n)

print 'Initial values [x,v]\n', X

ax = 0.0
ay = 0.0
az = 0.0

a = np.array([ax, ay, az])

nmax = int((tmax - t)/dt)

#-- define function to solve for x and v simultaneously
#-- ``acceleration a given by the grad of the potential function
#-- 10^6(2x^4 + 5y^4 + 10z^4)ms^2''; does that include the factor of
#-- q/m ?? Eqn 9 implies yes (a + v X B = 0)

def RHS(t, X):
	
#	output = np.zeros(n)
#	output[:,0] = X[:,1]
#	output[:,1] = 1.0
	output = np.zeros((h, n))
	output[:,0] = X[:,1]
	output[:,1] = a + (qe/me)*np.cross(X[:,1],b) 

	return output

rkX_arr = np.zeros((nmax, h, n))
#rkX_arr = np.zeros((nmax, n))

#-- rk4 solver

for ii in range(nmax):

	rkX_arr[ii,:,:] = X

	k1 = RHS(t, X)
	k2 = RHS(t + (dt/2.0), X + (dt/2.0)*k1)
	k3 = RHS(t + (dt/2.0), X + (dt/2.0)*k2)
	k4 = RHS(t + dt, X + dt*k3)

#	pdb.set_trace()

	X[:,1] += (dt/6.)*(k1[:,1] + 2.*k2[:,1] + 2.*k3[:,1] + k4[:,1])

	#-- define direction towards negative x (not really necessary, but done in the paper
	#-- that way)
	X[:,0] -= (dt/6.)*(k1[:,0] + 2.*k2[:,0] + 2.*k3[:,0] + k4[:,0])
	t += dt	

	#-- update acceleration with new position
	ax = 10.e6*(8.*X[0,0]**3.)
	ay = 10.e6*(20.*X[1,0]**3.)
	az = 10.e6*(40.*X[2,0]**3.)

	a = np.array([ax, ay, az])

t_arr = np.linspace(0, tmax-dt, nmax)

plt.figure(1)
plt.plot(rkX_arr[:,0,0], rkX_arr[:,0,1], label='vx')
plt.plot(rkX_arr[:,0,0], rkX_arr[:,1,1], label='vy')
plt.plot(rkX_arr[:,0,0], rkX_arr[:,2,1], label='vz')
plt.xlabel('Position (x)')
plt.ylabel('Velocity (ms$^{-1}$)')
plt.legend(loc='upper right', prop={'size':10})
plt.title('rk4')
plt.show(1)

#--------------------------------------#
# boris push method                    #
#--------------------------------------#

#-- electric field acceleration
#-- assuming the same as for rk4 (a contains all constant factors)
#-- otherwise there would be (qe/me) multiplying the aE expression as well. Easy enough to 
#-- add in if it's needed. 
aE = (dt / 2.0)*a

#-- rotation from B field
#-- (qB_mag/m)*(dt/2.)*e_para, b=B_mag*e_para
OM = (qe / me)*(dt / 2.0)*b
OM_mag_sq = OM[0]**2. + OM[1]**2. + OM[2]**2.
s = 2.*OM / (1.0 + OM_mag_sq)

boX_arr = np.zeros((nmax, h, n))

for ii in range(nmax):

	boX_arr[ii,:,0] = x
	boX_arr[ii,:,1] = v
	
	ax = 10.e6*(8.*x[0]**3.)
	ay = 10.e6*(20.*x[1]**3.)
	az = 10.e6*(40.*x[2]**3.)
       
	a = np.array([ax, ay, az])
	
	#-- update acceleration using new position values
	aE = (qe / me)*(dt / 2.0)*a
	
	vm = v + aE 
#	vp = vm + s*np.cross(v, OM) 
	vpr = vm + np.cross(vm, OM) 
	vpl = vm + np.cross(vpr, s) 
	v = vpl + aE 

	x = x - v*dt


#fig = plt.figure(2)
#ax = fig.gca(projection='3d')
#ax.plot(boX_arr[:,0,0], boX_arr[:,1,0], boX_arr[:,2,0])
#ax.set_title('Position (B0 = Bz, E0 = Ey, v0 = (vx,vy,vz))')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#plt.show(2)

plt.figure(3)
plt.plot(boX_arr[:,0,0], boX_arr[:,0,1], label='vx')
plt.plot(boX_arr[:,0,0], boX_arr[:,1,1], label='vy')
plt.plot(boX_arr[:,0,0], boX_arr[:,2,1], label='vz')
plt.xlabel('Position (x)')
plt.ylabel('Velocity (ms$^{-1}$)')
plt.legend(loc='upper right', prop={'size':11})
plt.title('Boris push')
plt.show(3)

v_perp1 = np.sum(R[0,:]*boX_arr[:,:,1], axis=1)
v_perp2 = np.sum(R[1,:]*boX_arr[:,:,1], axis=1)
v_para = np.sum(R[2,:]*boX_arr[:,:,1], axis=1)

plt.figure(4)
plt.plot(boX_arr[:,0,0], v_perp1, label='v$_{\perp,1}$')
plt.plot(boX_arr[:,0,0], v_perp2, label='V$_{\perp,2}$')
plt.plot(boX_arr[:,0,0], v_para, label='v$_{\parallel}$')
plt.xlabel('Position (x)')
plt.ylabel('Velocity (ms$^{-1}$)')
plt.legend(loc='upper right', prop={'size':11})
plt.title('Boris push')
plt.show(4)






























