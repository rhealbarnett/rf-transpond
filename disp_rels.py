#-----------------------------------------------------------#
# code CMA diagram                                          #
# x-axis om_p/om or density                                 #
# y-axis om_c/om or magnetic field                          #
# rlb 04012017                                              #
#-----------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import pdb

# --------------------------------------------------------- #

#-- setup variable density & magnetic field 

#r_max = 1.2
#r_min = 0.
#npts = 121
#r = np.linspace(r_min,r_max,npts+1)
#
##-- magnetic field of 1T
#B0_max = 3.4
#B0 = np.zeros(npts)
#for ii in range(0,npts):
#	B0[ii] = 1.0/(r_max - r[ii])
#	print ii, r[ii]
##B0[npts-1] = max(B0)
#B0 = (B0/max(B0))*B0_max
##pdb.set_trace()
#
#r = r[0:npts]

B0 = 1.0

#-- Swanson assumes mi/me = 5
#me = 9.11e-31
#mi = 2.*1.67e-27
me = 1.0
mi = 5.*me

#-- constants
#eps0 = 8.85e-12
#q_e = -1.602e-19
eps0 = 1.0
q_e = -1.0
q_i = abs(q_e)
sign_e = q_e / abs(q_e)
sign_i = q_i / abs(q_i)

#-- assume Ni=Ne=N
#-- from Swanson
#-- low density case om_pe^2/om_ce^2 = 0.3
#-- high density case om_pe^2/om_ce^2 = 3.0
#N_min = 1.0e16
#N_max = 1.0e18
#gradient = (N_min - N_max)/(r[0] - r[npts-1])
#N = gradient*(r - r[npts-1]) + N_max
#sig_sq = 0.2
#c = np.sqrt(sig_sq)
#b = 0.
#a = N_max/(c*np.sqrt(2.*np.pi))
#N = a*np.exp(((r - b)**2.)/(2.*c**2.))
#alpha = 35.
#N = np.exp(alpha*r) + N_min
alpha = 3.0
N = alpha*eps0*(B0**2.)/(me)
#N = 1.0e18

#-- plasma frequencies (assume ion has charge +e)
om_pe = np.sqrt((N*abs(q_e)**2)/(eps0*me))
om_pi = np.sqrt((N*abs(q_i)**2)/(eps0*mi))

#-- cycloctron frequencies
om_ce = abs(q_e)*B0/me
om_ci = abs(q_i)*B0/mi

#-- frequency values
fact = np.logspace(-3,1,100)
fact += 0.1*np.min(fact)
om = fact*om_ce
#freq = 51.0e6
#om = 2.*np.pi*freq

#-- solve dispersion relation n^2
#ns_R = 1.0 - ((om_pi**2)/(om*(om + om_ci))) - ((om_pe**2)/(om*(om - om_ce)))

S = 1.0 - ((om_pe**2.)/(om**2. - om_ce**2.) + (om_pi**2.)/(om**2. - om_ci**2.))
D = (sign_e*om_ce*om_pe**2.)/(om*(om**2. - om_ce**2.)) + (sign_i*om_ci*om_pi**2)/(om*(om**2. - om_ci**2.))
P = 1.0 - ((om_pe**2. / om**2.) + (om_pi**2. / om**2.))
#R = 1.0 - (((om_pe**2.)/(om*(om + sign_e*om_ce))) + ((om_pi**2.)/(om*(om + sign_i*om_ci))))
#L = 1.0 - (((om_pe**2.)/(om*(om - sign_e*om_ce))) + ((om_pi**2.)/(om*(om - sign_i*om_ci))))
R = S + D
L = S - D

S = S + 0.0j
D = D + 0.0j
P = P + 0.0j


dielec_tens = np.array([[S, -1.0j*D, 0.0],
			[1.0j*D, S, 0.0],
			[0.0, 0.0, P]])

theta = np.pi/2.0 

A = S*np.sin(theta)**2. + P*np.cos(theta)**2.
B = R*L*np.sin(theta)**2. + P*S*(1.0 + np.cos(theta)**2.)
C = P*R*L
#F = np.sqrt(((R*L - P*S)**2.)*np.sin(theta)**4. + 4.*P**2.*D**2.*np.cos(theta)**2.)
F = np.sqrt(B**2. - 4.*A*C)



ns_pos = (B + F) / (2.*A)
ns_neg = (B - F) / (2.*A)

pdb.set_trace()

n_pos_pos = np.sqrt(ns_pos)
n_pos_neg = -np.sqrt(ns_pos)
n_neg_pos = np.sqrt(ns_neg)
n_neg_neg = -np.sqrt(ns_neg)

plt.figure(1)
plt.plot(om/om_ce, np.real(n_pos_pos), label='Re(++)')
plt.plot(om/om_ce, np.real(n_pos_neg), label='Re(+-)')
plt.plot(om/om_ce, np.real(n_neg_pos), label='Re(-+)')
plt.plot(om/om_ce, np.real(n_neg_neg), label='Re(--)')
#plt.xlabel('r', size=16)
plt.axvline(om_pe, color='black', linestyle='--')
plt.axhline(1.0, color='black', linestyle='--')
plt.ylabel('$n$', size=16)
plt.xlabel('$\omega/\omega_{ce}$', size=16)
plt.xscale('log')
plt.ylim(-20,30)
plt.legend(loc='upper right')
plt.show(1)


plt.figure(2)
plt.plot(om/om_ce, np.imag(n_pos_pos), label='Im(++)')
plt.plot(om/om_ce, np.imag(n_pos_neg), label='Im(+-)')
plt.plot(om/om_ce, np.imag(n_neg_pos), label='Im(-+)')
plt.plot(om/om_ce, np.imag(n_neg_neg), label='Im(--)')
#plt.xlabel('r', size=16)
plt.axvline(om_pe, color='black', linestyle='--')
plt.axhline(1.0, color='black', linestyle='--')
plt.ylabel('$n$', size=16)
plt.xlabel('$\omega/\omega_{ce}$', size=16)
plt.xscale('log')
plt.ylim(-20,30)
plt.legend(loc='upper right')
plt.show(2)

plt.figure(3)
plt.plot(om/om_ce, ns_pos, label='+')
plt.plot(om/om_ce, ns_neg, label='-')
plt.ylabel('$n^2$')
plt.xlabel('$\omega/omega_{ce}$')
plt.axvline(om_pe, color='black', linestyle='--')
plt.axhline(1.0, color='black', linestyle='--')
plt.xscale('log')
plt.ylim(-20,30)
plt.legend(loc='upper right')
plt.show(3)
