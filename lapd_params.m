% ------------------------------------------------------------------- %
% LAPD parameter file for wave solver 
% 250219 rlbarnett, c3149416
% ------------------------------------------------------------------- %

% call constants
const = constants();

mp = const.mp;
me = const.me;
amu = const.amu;
e = -const.e;
q = const.e;
eps0 = const.eps0;
c0 = const.c0;

charge = [e; q];

% magnetic field magnitude, 1000 G = 0.1 T
B0 = 0.1;

% electron and ion temperatures, (eV)
% multiply by e for temps in K
Te = 1.0;
Ti = 6.0;

% ion mass : there are 3 possible ions in LAPD, 
% He, Ne and Ar. (kg)
mhe = 4.00*amu;
mne = 20.18*amu;
mar = 39.95*amu;

m = [me; mhe];

% plasma column is ~ 19 (m) 
xmin = -9.5;
xmax = 9.5;
npts = 512;

% driving frequency of the single strap, high power antenna (Hz)
freq = 2.4e6;
om = freq*2.0*pi;

% electron density range is (1.0e17 <= n <= 7.9e18) (m^-3)
Nmax = 7.9e18;
Nmin = 1.0e16;
n_new = logspace(log10(Nmin),log10(Nmax),npts);

% perpendicular wavenumber : just an approximation for now
% see figure 10 in Martin 2016 poster for n_perp and
% n_para, n = c0*k/om
n_perp = linspace(0,800,npts);
k_perp = om*n_perp./c0;
wave_perp = 2.0*pi./k_perp;
