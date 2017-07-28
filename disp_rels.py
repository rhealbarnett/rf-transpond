#-----------------------------------------------------------#
# code CMA diagram                                          #
# x-axis om_p/om or density                                 #
# y-axis om_c/om or magnetic field                          #
# rlb 04012017                                              #
#-----------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt

#%%

# --------------------------------------------------------- #
# om_p --> plasma frequency = ((n0*e^2)/(eps0*m))^(1/2)
# om_c --> electron cyclotron frequency == (|q|*B)/m
# v_th --> thermal velocity == sqrt((2*K*T_e)/m) 
# v_s  --> sound speed in plasma == (gamma_e*K*T_e + gamma_i*K*T_i)
# OM_c --> ion cyclotron frequency == (|q|*B)/M
# om_h --> upper hybrid frequency == om_p^2 + om_c^2
# va   --> Alfven velocity == B/sqrt(mu_0*rho)
# --------------------------------------------------------- #

# --------------------------------------------------------- #
# --- electromagnetic electron waves --- #
# -- light waves
# B0 = 0 : om^2 = om_p^2 + k^2*c^2
# -- O waves
# k _|_ B0, E1 || B0 : (c^2*k^2)/om^2 = 1 - om_p^2/om^2
# -- X wave
# k _|_ B0, E1 _|_ B0 : (c^2*k^2)/om^2 = 1 - (om_p^2/om^2)*((om^2 - om_p^2)/(om^2 - om_h^2))
# -- R wave (whistler mode)
# k || B0 : (c^2*k^2)/om^2 = 1 - (om_p^2/om^2)/(1 - (om_c^2/om^2))
# -- L wave
# k || B0 : (c^2*k^2)/om^2 = 1 - (om_p^2/om^2)/(1 + (om_c^2/om^2))
# --------------------------------------------------------- #

# --------------------------------------------------------- #
# --- electromagnetic ion waves --- #
# -- Alfven wave
# k || B0 : om^2 = k^2*va^2
# -- Magnetosonic wave
# k _|_ B0 : om^2/k^2 = c^2*((v_s^2 + va^2)/(c^2 + va^2))
# --------------------------------------------------------- #

import pdb

r_max = 1.2
r_min = 0.
npts = 121
r = np.linspace(r_min,r_max,npts+1)

#-- magnetic field of 1T
B0_max = 3.4
B0 = np.zeros(npts)
for ii in range(0,npts):
	B0[ii] = 1.0/(r_max - r[ii])
	print ii, r[ii]
#B0[npts-1] = max(B0)
B0 = (B0/max(B0))*B0_max
#pdb.set_trace()

r = r[0:npts]

#-- assume mi/me = 5
me = 9.11e-31
mi = 2.*1.67e-27
#freq = 51.0e6
#om = 2.*np.pi*freq

#-- constants
eps0 = 8.85e-12
q_e = -1.602e-19
q_i = abs(q_e)
sign_e = q_e / abs(q_e)
sign_i = q_i / abs(q_i)

#-- assume Ni=Ne=N
N_max = 1.0e18
N_min = 1.0e16
#gradient = (N_min - N_max)/(r[0] - r[npts-1])
#N = gradient*(r - r[npts-1]) + N_max
#sig_sq = 0.2
#c = np.sqrt(sig_sq)
#b = 0.
#a = N_max/(c*np.sqrt(2.*np.pi))
#N = a*np.exp(((r - b)**2.)/(2.*c**2.))
alpha = 35.
N = np.exp(alpha*r) + N_min
#alpha = 3.
#N = alpha*eps0*(B0**2.)/(me)

#-- plasma frequencies (assume ion has charge +e)
om_pe = np.sqrt((N*abs(q_e)**2)/(eps0*me))
om_pi = np.sqrt((N*abs(q_i)**2)/(eps0*mi))

#-- cycloctron frequencies
om_ce = abs(q_e)*B0/me
om_ci = abs(q_i)*B0/mi

#-- frequency values
#fact = np.logspace(-2,1, 100)
#fact += 0.1*np.min(fact)
#om = fact*om_ce
freq = 51.0e6
om = 2.*np.pi*freq

#-- solve dispersion relation n^2
ns_R = 1.0 - ((om_pi**2)/(om*(om + om_ci))) - ((om_pe**2)/(om*(om - om_ce)))

S = 1.0 - ((om_pe**2.)/(om**2. - om_ce**2.) + (om_pi**2.)/(om**2. - om_ci**2.))
D = (sign_e*om_ce*om_pe**2.)/(om*(om**2. - om_ce**2.)) + (sign_i*om_ci*om_pi**2)/(om*(om**2. - om_ci**2.))
P = 1.0 - ((om_pe**2. / om**2.) + (om_pi**2. / om**2.))
#R = 1.0 - (((om_pe**2.)/(om*(om + sign_e*om_ce))) + ((om_pi**2.)/(om*(om + sign_i*om_ci))))
#L = 1.0 - (((om_pe**2.)/(om*(om - sign_e*om_ce))) + ((om_pi**2.)/(om*(om - sign_i*om_ci))))
R = S + D
L = S - D

theta = np.pi/2.

A = S*np.sin(theta)**2. + P*np.cos(theta)**2.
B = R*L*np.sin(theta)**2. + P*S*(1.0 + np.cos(theta)**2.)
C = P*R*L
F = np.sqrt(((R*L - P*S)**2.)*np.sin(theta)**4. + 4.*P**2.*D**2.*np.cos(theta)**2.)
#F = np.sqrt(B**2. - 4.*A*C)

ns_pos = (B + F) / (2.*A)
ns_neg = (B - F) / (2.*A)
om_hybe_sq = (om_pe**2. + om_ce**2.)
om_hybi_sq = (om_pi**2. + om_ci**2.)
X_mode = ((om_pe**2)/(om**2.))*((om**2. - om_pe**2.)/(om**2. - om_hybe_sq))

plt.figure(1)
#plt.plot(om/om_ce, ns_R)
plt.plot(r, ns_pos, label='pos')
plt.plot(r, ns_neg, label='neg')
plt.xlabel('r', size=16)
plt.ylabel('$n^2$', size=16)
#plt.xscale('log')
plt.show(1)
