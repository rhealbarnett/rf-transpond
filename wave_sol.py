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

B0 = 3.4

#-- electron calcs
e = -1.6022e-19
me = 9.11e-31
om_ce = (e*B0) / me
Ne = 1.0e16
om_pe = np.sqrt((Ne*e**2.) / (eps0*me))

#-- ion calcs (assume fully ionised, 100% D (d) plasma, neutron mass approx proton mass)
q = 2.*abs(e)
mp = 1.67e-27
md = 2.*mp
om_cd = (q*B0) / md
Nd = 1.0e16
om_pd = np.sqrt((Nd*q**2.) / (eps0*me))

#-- wavenumbers (ky and kz from section IV in DVE 2015)
ky = 5.0
kz = 6.0

#-- non-rotated, cold plasma dielectric tensor 
S = 1.0 - om_pe**2/(om**2 - om_ce**2)
D = om_ce*om_pe/(om*(om**2 - om_ce**2))
P = 1.0 - om_pe**2/om**2
R = S + D
L = S - D

dielec_tens = np.array([[S, -1.0j*D, 0.0],
			[1.0j*D, S, 0.0],
			[0.0, 0.0, P]])

ns_perp = (S**2. + D**2.) / (S**2.)

#-- calculate kx for FW using dispersion relation
kx = np.sqrt((ns_perp*om**2./c0**2.) - ky**2. - kz**2. + 0.0j)

#-- calculate kx for SW using dispersion relation


#-- rotation matrix
alpha = 0.01
beta = 0.01

rot = np.array([[np.cos(beta)*np.cos(alpha), np.cos(beta)*np.sin(alpha), -np.sin(beta)],
		[-np.sin(alpha), np.cos(alpha), 0.0],
		[np.sin(beta)*np.cos(alpha), np.sin(beta)*np.sin(alpha), np.cos(beta)]])

#-- dispersion relation matrix
disp_rel = np.array([[ky**2. + kz**2. - k0**2.*S, -ky*kx + k0**2.*1.0j*D, -kz*kx],
		     [-ky*kx - k0**2.*1.0j*D, kz**2. + kx**2. - k0**2.*S, -kz*ky],
		     [-kz*kx, -ky*kz, ky**2. + kx**2. - k0**2.*P]])	

eigvals, eigvecs = np.linalg.eig(disp_rel)

eig_check = np.array([[1.0, 1.0j],
		      [-1.0j, 1.0]])

vals, vecs = np.linalg.eig(eig_check)

