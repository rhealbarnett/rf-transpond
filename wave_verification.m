% ------------------------------------------------------------------- %
% Wave verification parameter file for wave solver 
% 100819 rlbarnett, c3149416
% ------------------------------------------------------------------- %

%--
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

%--
% Magnetic field magnitude.
B0 = 2.5;

%-- 
% Set some domain size and length
xmin = -4.;
xmax = 4.0;
npts = 512;
xax = linspace(xmin,xmax,npts);

%--
% Particle masses
me = me*ones(1,npts);
mp = mp*ones(1,npts);

% DLG - since I don't have the license for "makedist"
% I fixed your cos ramping function :)
% dampFac = 25.0e1;
% np_bound = floor(0.2*npts);
% ax = linspace(0,pi,np_bound);
% damp0 = (cos(ax)+1)/2;
% damp = ones(1,npts);
% damp(1:np_bound) = damp(1:np_bound) + dampFac*1i*damp0;
% damp(end-np_bound+1:end) = damp(end-np_bound+1:end) + dampFac*1i*fliplr(damp0);
% 
% me = me .* damp;

m_s = [me; mp];

%--
% Driving frequency in the MHz range
freq = 51.0e6;
om = freq*2.0*pi;

%--
% Scan over density values.
Nmax = 1.0e20;
Nmin = 1.0e15;
% n_new = logspace(log10(Nmin),log10(Nmax),npts);
n_new = 1.0e19*ones(1,npts);

%--
% Calculate free space wavenumber and set perpendicular (or parallel) 
% wave number. 
k0 = om/c0;
k_para = 5.0;
ky = 0.;
kz = k_para;

%--
% Set current source parameters. 
source_width = 0.06;
source_loc = 0;
source_mult = 1.0;

