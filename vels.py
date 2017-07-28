#---------------------------------------------#
# electron velocities 			      #
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
import math

T_ev = 5.0 
B_mag = 3.4 
#m = 9.11e-31
u = 1.66053940e-27
m = 3.0160293*u
e = 1.6022e-19
kb = 1.38e-23
vt = np.sqrt((T_ev*e) / m)
q = 2.*e
om_c = (q*B_mag) / m
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

e = np.array([R[0,:], R[1,:], R[2,:]])
#e = np.array([0.0,0.0,1.0])
e_para = e[2]

b = e_para*B_mag

#-- initial velocity
#- electron initial velocity
v0 = vt/3.
#v0 = 0.1*vt

vx = v0
vy = v0
vz = v0
v = np.array([vx, vy, vz])

rl = m*v[1]/(q*b[2])

#-- initial position
xx = 0.0
xy = 0.0
xz = 0.0
x = np.array([xx, xy, xz])

t = 0.0
num_cyc = 10
num_points = 100
dt = period/num_points
tmax = num_cyc*num_points*dt
nmax = int((tmax - t)/dt)

#-- position and velocity
n = 2
#-- number of spatial (position) dimensions
pos_dim = 3

#---------------------------------------#
# RK4    			        #
#---------------------------------------#

X = np.zeros((pos_dim, n))
X[:,0] = x
X[:,1] = v

print 'Initial values [x,v]\n', X

ax = 0.0
ay = 0.0
az = 0.0

a = np.array([ax, ay, az])

#-- define function to solve for x and v simultaneously
#-- ``acceleration a given by the grad of the potential function
#-- 10^6(2x^4 + 5y^4 + 10z^4)ms^2''; does that include the factor of
#-- q/m ?? Eqn 9 implies yes (a + v X B = 0)

#def RHS(t, X):
#	
##	output = np.zeros(n)
##	output[:,0] = X[:,1]
##	output[:,1] = 1.0
#	output = np.zeros((h, n))
#	output[:,0] = X[:,1]
#	output[:,1] = a + (qe/m)*np.cross(X[:,1],b) 
#
#	return output
#
#rkX_arr = np.zeros((nmax, h, n))
##rkX_arr = np.zeros((nmax, n))
#
##-- rk4 solver
#
#for ii in range(nmax):
#
#	rkX_arr[ii,:,:] = X
#
#	k1 = RHS(t, X)
#	k2 = RHS(t + (dt/2.0), X + (dt/2.0)*k1)
#	k3 = RHS(t + (dt/2.0), X + (dt/2.0)*k2)
#	k4 = RHS(t + dt, X + dt*k3)
#
##	pdb.set_trace()
#
#	X[:,1] += (dt/6.)*(k1[:,1] + 2.*k2[:,1] + 2.*k3[:,1] + k4[:,1])
#
#	#-- define direction towards negative x (not really necessary, but done in the paper
#	#-- that way)
#	X[:,0] -= (dt/6.)*(k1[:,0] + 2.*k2[:,0] + 2.*k3[:,0] + k4[:,0])
#	t += dt	
#
#	#-- update acceleration with new position
#	ax = 10.e6*(8.*X[0,0]**3.)
#	ay = 10.e6*(20.*X[1,0]**3.)
#	az = 10.e6*(40.*X[2,0]**3.)
#
#	a = np.array([ax, ay, az])
#
#t_arr = np.linspace(0, tmax-dt, nmax)
#
#plt.figure(1)
#plt.plot(rkX_arr[:,0,0], rkX_arr[:,0,1], label='vx')
#plt.plot(rkX_arr[:,0,0], rkX_arr[:,1,1], label='vy')
#plt.plot(rkX_arr[:,0,0], rkX_arr[:,2,1], label='vz')
#plt.xlabel('Position (x)')
#plt.ylabel('Velocity (ms$^{-1}$)')
#plt.legend(loc='upper right', prop={'size':10})
#plt.title('rk4')
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.show(1)

#--------------------------------------#
# boris push method                    #
#--------------------------------------#

#-- electric field acceleration
#-- assuming the same as for rk4 (a contains all constant factors)
#-- otherwise there would be (qe/me) multiplying the aE expression as well. Easy enough to 
#-- add in if it's needed.

pos_prev = np.array([0.0,0.0,0.0])

#-- rotation from B field
#-- (qB_mag/m)*(dt/2.)*e_para, b=B_mag*e_para
OM = (q / m)*(dt / 2.0)*b
OM_mag_sq = OM[0]**2. + OM[1]**2. + OM[2]**2.
s = 2.*OM / (1.0 + OM_mag_sq)

boX_arr = np.zeros((nmax, pos_dim, n))
a_arr = np.zeros((nmax, pos_dim))
v_drift = np.zeros((num_cyc, pos_dim))

count = 0
pos_av = np.zeros(pos_dim)
pos_avarr = np.zeros((num_cyc, pos_dim))
a_avarr = np.zeros((num_cyc, pos_dim))
v_mag_arr = np.zeros((nmax, pos_dim))
print v

for ii in range(nmax):

	v_mag = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
	v_mag_arr[ii,:] = v_mag
	vp1 = np.dot(v,R[0,:])
	vp2 = np.dot(v,R[1,:])
	vpara = np.dot(v,e_para)
	v_check = np.sqrt(vp1**2 + vp2**2 + vpara**2)
#	print 'v magnitude', v_mag
#	print 'vp1, vp2, vpara', vp1, vp2, vpara
#	print 'v check', v_check
	

	boX_arr[ii,:,0] = x
	boX_arr[ii,:,1] = v
	
	#-- update acceleration using new positions	
	ax = 1.e6*(8.*x[0]**3.)
	ay = 1.e6*(20.*x[1]**3.)
	az = 1.e6*(40.*x[2]**3.)
 
	a = np.array([ax, ay, az])
	a_arr[ii,:] = a	

	aE = (dt / 2.0)*a
	
	vm = v + aE 
	vpr = vm + np.cross(vm, OM)
	vpl = vm + np.cross(vpr, s)
	v = vpl + aE 

	x = x + v*dt
	pos_av = pos_av + x

	#-- 'brute force' drift velocity calculation
	#-- average position over a cyclotron period
	if math.fmod(ii,num_points)==0.0:
		xa = pos_av/num_points
		pos_avarr[count,:] = xa
		ds = np.sqrt((xa - pos_prev)**2)
		d_vel = ds/period
	#	print "E", a
	#	print "B", b
	#	print "v drift", d_vel
		v_drift[count,:] = d_vel
		a_avarr[count,:] = a
		pos_prev = xa
		pos_av = np.zeros(pos_dim)		
		count = count + 1


#-- rotate velocities to perp and parallel components
v_perp1 = np.sum(R[0,:]*boX_arr[:,:,1], axis=1)
v_perp2 = np.sum(R[1,:]*boX_arr[:,:,1], axis=1)
v_para = np.sum(R[2,:]*boX_arr[:,:,1], axis=1)
vd_perp1 = np.sum(R[0,:]*v_drift, axis=1)
vd_perp2 = np.sum(R[1,:]*v_drift, axis=1)
vd_para = np.sum(R[2,:]*v_drift, axis=1)

#-- calculate dv_(drift)/dt  
dv_dt = (v_drift[1:num_cyc,:] - v_drift[0:num_cyc-1,:])/period

#-- v drift from eqn 8 from DVE 2015 paper, page 5
vd = (1.0/om_c)*(np.cross(a_avarr, e_para))
vcom = np.cross(vd, om_c*e_para)

vd_mag = v_mag_arr*e_para
vd_mag_perp1 = np.sum(R[0,:]*vd_mag, axis=1)
vd_mag_perp2 = np.sum(R[1,:]*vd_mag, axis=1)
vd_mag_para = np.sum(R[2,:]*vd_mag, axis=1)

plt.figure(2)
plt.plot(boX_arr[:,0,0], boX_arr[:,0,1], label='vx')
plt.plot(boX_arr[:,0,0], boX_arr[:,1,1], label='vy')
plt.plot(boX_arr[:,0,0], boX_arr[:,2,1], label='vz')
plt.plot(boX_arr[:,0,0], vd_mag[:,0], '--', label='v$_{d,x}$')
plt.plot(boX_arr[:,0,0], vd_mag[:,1], '--',  label='v$_{d,y}$')
plt.plot(boX_arr[:,0,0], vd_mag[:,2], '--', label='v$_{d,z}$')
plt.xlabel('Position ($m$)', fontsize=16)
plt.ylabel('Velocity ($ms^{-1}$)', fontsize=16)
plt.legend(loc='lower right', prop={'size':10})
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.title('e$^-$, T=%r$eV$' % (T_ev))
plt.show(2)

plt.figure(3)
plt.plot(boX_arr[:,0,0], v_perp1, label='v$_{\perp,1}$')
plt.plot(boX_arr[:,0,0], v_perp2, label='v$_{\perp,2}$')
plt.plot(boX_arr[:,0,0], v_para, label='v$_{\parallel}$')
plt.plot(boX_arr[:,0,0], vd_mag_perp1, '--', label='v$_{\perp,1,d}$')
plt.plot(boX_arr[:,0,0], vd_mag_perp2, '--',  label='v$_{\perp,2,d}$')
plt.plot(boX_arr[:,0,0], vd_mag_para, '--', label='v$_{\parallel,d}$')
plt.xlabel('Position ($x$)')
plt.ylabel('Velocity ($ms^{-1}$)')
plt.legend(loc='upper left', prop={'size':10})
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.title('e$^-$, T=%r$eV$' % (T_ev))
plt.show(3)

plt.figure(4)
plt.plot(pos_avarr[:,0], vcom[:,0], 'r--',  label='[v$_s\\times\Omega$]$_x$')
plt.plot(pos_avarr[:,0], a_avarr[:,0],'g--', label='a$_x$')
plt.plot(pos_avarr[:,0], a_avarr[:,0]+vcom[:,0], 'c--', label='[a + v$_s\\times\Omega$]$_x$')
plt.xlabel('Position ($m$)', fontsize=16)
plt.ylabel('Acceleration ($ms^{-2}$)', fontsize=16)
plt.legend(loc='upper left')
plt.title('$^3$He, T=%r$eV$' % (T_ev))
plt.xlim(min(pos_avarr[:,0]), max(pos_avarr[:,0]))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show(4)










