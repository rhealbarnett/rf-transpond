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
npts = 16;

% driving frequency of the single strap, high power antenna (Hz)
freq = 2.4e6;
om = freq*2.0*pi;

% density range is (1.0e10 <= n <= 2.0e12) (cm^-3)
% which is (1.0e16 <= n <= 2.0e18) (m^-3)
Nmax = 1.0e17;
n = Nmax*ones(1,npts);


