#-----------------------------------#
# Solve simple Poissons equation    #
# rlb 150617 			    #
#-----------------------------------#

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.special import erf
import pdb


npts = 101

xmin = 0
xmax = 0.2
ymin = 0
ymax = 100

dx = (xmax - xmin)/(npts - 1)
dy = (ymax - ymin)/(npts - 1)

x = np.linspace(xmin, xmax, npts)
y = np.linspace(ymin, ymax, npts)

Nemax = 5.0
Ne = np.exp(203.8*x)
Nimax = 2.0
Ni = np.exp(199.2*x)
e = 1.6022e-19
eps0 = 8.85e-12

#rho = -(e / eps0)*(Ne - Ni)
rho = np.ones(npts)*1.0e4
#rho = np.zeros(npts)

#phi = np.zeros((npts, npts))
#phi[0,:] = 0.0
#phi[npts-1,:] = 0.0
#phi[:,0] = 0.0
#phi[:,npts-1] = 1.0

#--- poisson

phi = np.zeros(npts)
phi[0] = 0.
phi[npts-1] = 0.
#pdb.set_trace()
phi_old = np.zeros(npts)

threshold = 1.0e-5
opt = 2.0 / (1.0 + (np.pi/npts))
nmax = 100000 
#rms_err = np.zeros(nmax)

for ii in range(nmax):

#	phi[1:npts-1,1:npts-1] = (1.0/4.0)*(phi[0:npts-2,1:npts-1] + phi[2:npts,1:npts-1] + phi[1:npts-1,0:npts-2] + phi[1:npts-1,2:npts]) 
	phi[1:npts-1] = (1.0/2.0)*(phi[0:npts-2] + phi[2:npts] - (dx**2.)*rho[1:npts-1])

#	pdb.set_trace()

	rms_err = np.sqrt(np.sum((phi - phi_old)**2)/npts)

	phi_old[:] = phi

	if rms_err<=threshold:
		break
	else:
		continue	
#
#	pdb.set_trace()

#fig = plt.figure(1)
#ax = fig.gca(projection='3d')
#X, Y = np.meshgrid(x, y)
#surf = ax.plot_surface(X, Y, phi, cmap=cm.coolwarm, linewidth=0)
#ax.set_zlim(0., 1.)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.show(1)

#--- wave eqn

alpha = 0.5
beta = 0.5



R = np.array([[np.cos(beta)*np.cos(alpha), np.cos(beta)*np.sin(alpha), -np.sin(beta)],
		[-np.sin(alpha), np.cos(alpha), 0],
		[np.sin(beta)*np.cos(alpha), np.sin(beta)*np.sin(alpha), np.cos(beta)]])




