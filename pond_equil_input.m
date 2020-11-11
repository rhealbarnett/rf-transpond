%------------------------------------------------------------------------%
% pond_equil_input.m 
%
% Input parameter script for ponderomotive equilibrium density 
% calculation.
%
% rlbarnett 201111
%------------------------------------------------------------------------%

% ------------------------------------------------------------------- %
% LAPD parameter file for wave solver 
% 250219 rlbarnett, c3149416
% Updated to reflect LAPD RF experimental campaign data 050919
% Updated for cpc manuscript 201021
% ------------------------------------------------------------------- %

%--
% Import physical constants.
const = constants();

mp = const.mp;
me = const.me;
amu = const.amu;
e = -const.e;
q = const.e;
eps0 = const.eps0;
c0 = const.c0;
mu0 = const.mu0;

%--
% Actual plasma column is ~ 18 (m). However, use reduced size as interest
% is close to the antenna. 
zmin = -4.0;
zmax = 4.0;

%-- 
% npts will be set with the transport equilibrium script.
npts = 512;
zax = linspace(zmin,zmax,npts);
dz = (zmax - zmin) / (npts);

%--
% Driving frequency of the single strap, high power antenna (Hz)
% Driven at 2.38MHz, but FFT of experimental data shows it is closer to
% ~2.52MHz. 
freq = 2.52e6;
om = freq*2.0*pi;

%--
% Ion mass : He only
mhe = 4.00*amu;
mhe = mhe*ones(1,npts);

%--
% Electron mass.
me = me*ones(1,npts);

%--
% Collect masses in 2d array for dielectric tensor calculation. 
m_s = [me; mhe];

%--
% Thermal velocity 
Te = 5.0;
Ti = 0.5;
v_th = sqrt((Te + Ti)*abs(e)/mhe(1));

%--
% Initial density, flat.
n_init = 1.0e17*ones(1,npts);

%--
% Spatial electric field profile (gaussian)
E_width = (0.5)/(2.*sqrt(2.*log(2.)));
E_loc = 0.0;
E_mult = 1.0e9;

gauss_mult = 1.0/(E_width*sqrt(2.0*pi));
E_0 = gauss_mult*exp(-(zax - E_loc).^2/(2.0*E_width^2));
E_0 = E_0 / max(E_0);
E_0 = E_0*E_mult;





