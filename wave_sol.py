#------------------------------------#
# wave solver  			     #
# del X del X E = k0^2 K.E           #
# rlb 290617  			     #
#------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import pdb

freq = 51.e6
om = 2.*np.pi*freq
mu0 = 4.*np.pi*1.0e-7
eps0 = 8.85e-12
c0 = 1.0 / np.sqrt(mu0*eps0)
k0 = om / c0

B0 = 1.0
e = 1.6022e-19
me = 9.11e-31
om_ce = (e*B0) / me
Ne = 1.0e16
om_pe = np.sqrt((Ne*e**2.) / (eps0*me))

alpha = 0.1
beta = 0.

R = np.array([[np.cos(beta)*np.cos(alpha), np.cos(beta)*np.sin(alpha), -np.sin(beta)],
		[-np.sin(alpha), np.cos(alpha), 0.0],
		[np.sin(beta)*np.cos(alpha), np.sin(beta)*np.sin(alpha), np.cos(beta)]])




