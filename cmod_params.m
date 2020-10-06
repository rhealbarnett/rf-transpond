% ------------------------------------------------------------------- %
% Alcator C-Mod parameter file 
% 200728 rlbarnett, c3149416
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
% Magnetic field magnitude, 1000 G = 0.1 T. Assume it is fully aligned with 
% the z coordinate. 
B0 = 5.5;

%--
% Actual plasma column is ~ 18 (m). However, use reduced size as interest
% is close to the antenna. 
zmin = -4.0;
zmax = 4.0;
npts = 512;
zax = linspace(zmin,zmax,npts);
dz = (zmax - zmin) / (npts-1);

%--
% Driving frequency of the single strap, high power antenna (Hz)
% Driven at 2.38MHz, but FFT of experimental data shows it is closer to
% ~2.52MHz. 
freq = 75.e6;
om = freq*2.0*pi;
% om_i = 0.1i*om;

%--
% Ion mass : He only
md = 2.014*amu;
md = md*ones(1,npts);

%--
% Electron mass.
me = me*ones(1,npts);

%--
% Collect masses in 2d array for dielectric tensor calculation. 
m_s = [me; md];

%--
% Electron density range is (1.0e17 <= n <= 7.9e18) (m^-3). Scan over these
% values, +- some amount. 
% Nmin = 1.0e16;
% Nmax = 1.0e19;
% n_new = logspace(log10(Nmin),log10(Nmax),npts);
n_new = 1.e16*ones(1,npts);

%--
% Wavenumber in x approximated using experimental data, kx ~ (0 + 20i)
% m^-1. 
k0 = (om/c0);
kx = 165;
% kz = linspace(3,14,100);
ky = 0.0;
% ky = linspace(0,600,300);
k_perp = sqrt(kx.^2 + ky.^2); 
n_refrac = c0*k_perp./om;

%--
% Wavenumbers as function of spatial location (currently constant).
dampk = ones(1,npts);
kx = kx.*dampk;
ky = ky.*dampk;

%-- 
% Current source parameters.
source_width = (0.08)/(2.*sqrt(2.*log(2.)));
source_loc = 0.093;
source_mult = 9.065e6;
damp_len = 0.35;

mult = 1.0/(source_width*sqrt(2.0*pi));
source_right = mult*exp(-(zax - source_loc).^2/(2.0*source_width^2));

source_loc = -0.093;
mult = -1.0/(source_width*sqrt(2.0*pi));
source_left = mult*exp(-(zax - source_loc).^2/(2.0*source_width^2));

source_loc = 0.093 + 0.186;
mult = -1.0/(source_width*sqrt(2.0*pi));
source_r2 = mult*exp(-(zax - source_loc).^2/(2.0*source_width^2));

source_loc = -0.093 - 0.186;
mult = 1.0/(source_width*sqrt(2.0*pi));
source_l2 = mult*exp(-(zax - source_loc).^2/(2.0*source_width^2));

source = source_right + source_left + source_r2 + source_l2;
source = source / max(source);
source = source*source_mult;

