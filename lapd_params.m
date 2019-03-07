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
mu0 = const.mu0;

charge = [e; q];

% magnetic field magnitude, 1000 G = 0.1 T
B0 = 0.1;

% electron and ion temperatures, (eV)
% multiply by e for temps in K
Te = 5.0;
Ti = 0.5;

% plasma column is ~ 18 (m) 
xmin = -4.;
xmax = 4.0;
% npts = 2048;
xax = linspace(xmin,xmax,npts);

% ion mass : there are 3 possible ions in LAPD, 
% He, Ne and Ar. (kg)
mhe = 4.00*amu;
mne = 20.18*amu;
mar = 39.95*amu;

%--
% electron mass
me = me*ones(1,npts);
mhe = mhe*ones(1,npts);

% DLG - since I don't have the license for "makedist"
% I fixed your cos ramping function :)
dampFac = 50.0e1;
np_bound = floor(0.2*npts);
ax = linspace(0,pi,np_bound);
damp0 = (cos(ax)+1)/2;
damp = ones(1,npts);
damp(1:np_bound) = damp(1:np_bound) + dampFac*1i*damp0;
damp(end-np_bound+1:end) = damp(end-np_bound+1:end) + dampFac*1i*fliplr(damp0);

me = me .* damp;

m_s = [me; mhe];

% driving frequency of the single strap, high power antenna (Hz)
freq = 2.4e6;
% freq = 13.0e6;
om = freq*2.0*pi;

% electron density range is (1.0e17 <= n <= 7.9e18) (m^-3)
% Nmax = 7.9e18;
% Nmin = 1.0e17;
% n_new = logspace(log10(Nmin),log10(Nmax),npts);
% n_new = Nmin*ones(1,npts);

% perpendicular wavenumber : just an approximation for now
% see figure 10 in Martin 2016 poster for n_perp and
% n_para, n = c0*k/om
% n_perp = linspace(0,800,800);
% n_perp = 200;
% k_perp = om*n_perp./c0;
k_perp = 5.;
n_perp = c0*k_perp./om;
% wave_perp = 2.0*pi./k_perp;
% k0 = om/c0;
% wave0 = 2.0*pi/k0;

% kx = 2.0*pi/wave;
% kx = 40.;

ky = 0.;
kz = k_perp;
k0 = om/c0;

source_width = 0.03;
source_loc = 0;


