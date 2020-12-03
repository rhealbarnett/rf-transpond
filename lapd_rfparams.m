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
% Collect charge in a 2d array for dielectric tensor calculation.
q_s = [e; q];

%--
% LAPD magnetic field magnitude, 1000 G = 0.1 T. Assume it is aligned with 
% the z coordinate. 
B0 = 0.1;

%--
% Actual plasma column is ~ 18 (m). However, use reduced size as interest
% is close to the antenna. 
zmin = -4.0;
zmax = 4.0;

%-- 
% npts will be set with the transport equilibrium script.
zax = linspace(zmin,zmax,npts);
dz = (zmax - zmin) / (npts-1);

%--
% Driving frequency of the single strap, high power antenna (Hz)
% Driven at 2.38MHz, but FFT of experimental data shows it is closer to
% ~2.52MHz. 
freq = 2.38e6;
om = freq*2.0*pi;
period = 1.0/freq;

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
% Wavenumber in x approximated using experimental data, kx ~ (0 + 20i)
% m^-1. 
k0 = (om/c0);
kx = 20.0i;
ky = 0.0;
k_perp = sqrt(kx.^2 + ky.^2); 
n_refrac = c0*k_perp./om;

%--
% Wavenumbers as function of spatial location (currently constant).
dampk = ones(1,npts);
kx = kx.*dampk;
ky = ky.*dampk;

%-- 
% Current source parameters.
source_width = (0.06)/(2.*sqrt(2.*log(2.)));
source_loc = 0;

%--
% Scaling for the current source calculated from experimental data.
source_mult = 2.0e6;

%--
% Calculate Gaussian source term.
mult = 1.0/(source_width*sqrt(2.0*pi));
source = mult*exp(-(zax - source_loc).^2/(2.0*source_width^2));
source = source / max(source);
source = source*source_mult;

%--
% Number of points to include in absorbing boundary region.
damp_len = 0.2;
dampFac = 1.0e3;

%--
% Set transport switches
MMS = 0;
momentum = 0;
continuity = 0;
sfile = 0;
couple = 1;
test = 0;






