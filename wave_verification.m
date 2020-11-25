% ------------------------------------------------------------------- %
% Wave verification parameter file for wave solver 
% 100819 rlbarnett, c3149416
% ------------------------------------------------------------------- %

%--
% call constants
const = constants();

mp = const.mp;
me = const.me;
e = -const.e;
q = const.e;
eps0 = const.eps0;
c0 = const.c0;
mu0 = const.mu0;
amu = const.amu;

q_s = [e; q];

%--
% Magnetic field magnitude.
B0 = 2.5;

%--
% Driving frequency in the MHz range
freq = 52.0e6;
om = freq*2.0*pi;

%-- 
% Set some domain size and length
zmin = -4.;
zmax = 4.0;
npts = 1024;
zax = linspace(zmin,zmax,npts);
dz = (zmax - zmin) / (npts-1);

%--
% Particle masses
me = me*ones(1,npts);
mp = mp*ones(1,npts);
mhe = 4.0*amu;
mhe = mhe*ones(1,npts);

% m_s = [me; mp].*damp;
m_s = [me; mhe];

%--
% Scan over density values.
Nmax = 1.0e20;
Nmin = 1.0e19;
n_new = logspace(log10(Nmin),log10(Nmax),npts);

%--
% Calculate free space wavenumber and set perpendicular (or parallel) 
% wave number. 
k0 = om/c0;
kx = 20.;
ky = 0.;

k_perp = sqrt(kx.^2 + ky.^2); 
n_refrac = c0*k_perp./om;

dampk = ones(1,npts);
kx = kx*dampk;
ky = ky*dampk;


%--
% Set current source parameters. 
source_width = (0.06)/(2.*sqrt(2.*log(2.)));
source_loc = 0;
source_mult = 1.0;

mult = 1.0/(source_width*sqrt(2.0*pi));
source = mult*exp(-(zax - source_loc).^2/(2.0*source_width^2));
source = source / max(source);
source = source*source_mult;

damp_len = 0.3;


