% ------------------------------------------------------------------- %
% LAPD parameter file for wave solver 
% 250219 rlbarnett, c3149416
% Updated to reflect LAPD RF experimental campaign data 050919
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
B0 = 0.1;

%--
% Actual plasma column is ~ 18 (m). However, use reduced size as interest
% is close to the antenna. 
xmin = -4.;
xmax = 4.0;
% npts = 2048;
xax = linspace(xmin,xmax,npts);
dx = (xmax - xmin) / (npts-1);

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
% Electron density range is (1.0e17 <= n <= 7.9e18) (m^-3). Scan over these
% values, +- some amount. 
% Nmax = 1.0e19;
% Nmin = 1.0e16;
% n_new = logspace(log10(Nmin),log10(Nmax),npts);
% n_new = 1.0e17*ones(1,npts);

%--
% Wavenumber in x approximated using experimental data, kx ~ (0 + 20i)
% m^-1. 
k0 = (om/c0);
kx = 20.0i;
ky = 0.0;
% ky = linspace(0,40,100);
k_perp = sqrt(kx.^2 + ky.^2); 
n_perp = c0*k_perp./om;

%--
% Wavenumbers as function of spatial location (currently constant).
dampk = ones(1,npts);
kx = kx.*dampk;
ky = ky.*dampk;

%-- 
% Current source parameters.
source_width = 0.06/(2.*sqrt(2.*log(2.)));
source_loc = 0;
source_mult = 37000;

mult = 1.0/sqrt(2.0*pi*source_width);
source = mult*exp(-(xax - source_loc).^2/(2.0*source_width^2));
source = source / max(source);
source = source*source_mult;



