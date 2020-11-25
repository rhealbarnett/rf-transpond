%------------------------------------------------------------------------%
% pond_equil_input.m 
%
% Input parameter script for ponderomotive equilibrium density 
% calculation.
%
% rlbarnett 201111
%------------------------------------------------------------------------%

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
kb = const.kb;

%--
% Load transport equilibrium file
filepath = '/Volumes/DATA/LAPD/matlab/inputs/';
load(strcat(filepath,'equil_transport_input.mat'),'npts','nxax','vxax',...
    'cs','ndx','vdx','n_new','vx_new','n_source');

%--
% Set domain limits.
zmin = -4.0;
zmax = 4.0;

%-- 
% npts will be set with the transport equilibrium script.
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

q_s = [e; q];

%--
% Thermal velocity 
Te = 5.0;
Ti = 0.5;
v_th = sqrt((Te + Ti)*abs(e)/(me(1)));

%--
% Initial density.
n_init = n_new;
n_init_uni = interp1(nxax,n_init,zax,'linear');

%--
% Spatial electric field profile (gaussian)
E_width = (0.5)/(2.*sqrt(2.*log(2.)));
E_loc = 0.0;
E_mult = 2.0e2;

gauss_mult = 1.0/(E_width*sqrt(2.0*pi));
E_0 = gauss_mult*exp(-(zax - E_loc).^2/(2.0*E_width^2));
E_0 = E_0 / max(E_0);
E_0 = E_0*E_mult;

%--
% Calculate equilibrium density profile for given E field
n_final = pond_equil(E_0,n_init_uni,om,v_th,abs(e),(me(1)));

%--
% Assign transport parameters from equilibrium transport file to variables.
vx_init = vx_new;
LuBC = -cs;
RuBC = cs;

%--
% 
dt = 0.99*min(ndx)/cs;
m = mhe(1);
tmax = 2.0e-3;
nmax = round(tmax/dt);
nu = 1.0;
save_iter = round(nmax/10);

%--
%
rf_ex = zeros(1,npts-1);
rf_ey = zeros(1,npts-1);
rf_ez = E_0;
