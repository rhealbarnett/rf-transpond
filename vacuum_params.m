% ------------------------------------------------------------------- %
% LAPD parameter file for wave solver 
% 250219 rlbarnett, c3149416
% ------------------------------------------------------------------- %

% call constants
const = constants();

mp = const.mp;
me = const.me;
e = -const.e;
eps0 = const.eps0;
c0 = const.c0;
mu0 = const.mu0;

% magnetic field magnitude, 1000 G = 0.1 T
B0 = 0.0;

% plasma column is ~ 18 (m) 
xmin = 0.;
xmax = 1.;
% npts = 256;
xax = linspace(xmin,xmax,npts);

%--
% electron mass
me = me*ones(1,npts);

% zero density
n = zeros(1,npts);

% perpendicular wavenumber : none
Ex = 1.0;
Ey = 1.0;
Ez = 0.01;

nwave = 5.;
wave = (xmax - xmin)/(nwave);
kx = 2.0*pi/wave;
freq = 13.0e6;
om = 2.0*pi*freq;
ky = 0.;
kz = 0.;
k0 = om/c0;

ex_solx = Ex*exp(1i*kx*xax);
ex_soly = Ey*exp(1i*kx*xax);
ex_solz = Ez*exp(1i*kx*xax);

ex_sol = zeros(1,3*npts);
ex_sol(1:3:3*npts) = ex_solx;
ex_sol(2:3:3*npts) = ex_soly;
ex_sol(3:3:3*npts) = ex_solz;

source_width = 0.01;
        
        
