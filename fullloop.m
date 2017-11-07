%----------------------------------------%
% main script for 1D DVE loop            %
% rlbarnett c3149416 061117              %
%----------------------------------------%

clear all; close all

%--
% constants
mu0 = 4.0*pi*1.0e-7;
eps0 = 8.85e-12;
c0 = 1.0/sqrt(eps0*mu0);

%--
% magnetic field (tesla)
B0 = 2.5;

%--
% driver freq
freq = 51.0e6;
om = 2.0*pi*freq;
k0 = om/c0;
wavel0 = (2*pi)/k0;

%--
% define wavenumbers ky and kz (/m); use values given in van eester section IV?? No
% others mentioned
ky = 5.0;
kz = 6.0;

%--
% "common local derivatives for N0, v||^2, static potential and
% ponderomotive potential" 
lamby = 0.0;
lambz = 0.0;

%--
% spatial domain
npts = 1000;
dx = 0.002;
xmin = 0.0;
xmax = 0.2;
% npts = ((xmax - xmin)/dx);
xax = linspace(xmin, xmax, npts);
% xax = xmin:dx:xmax;
% dx = (xmax - xmin)/(npts - 1);

%--
% density -- set to zero for vacuum case
Nmax = 5.0e17;
N0 = Nmax*ones(npts,1);

%--
% electron constants
e = -1.6022e-19;
me = 9.11e-31;

%--
% temperature
T_ev = 15.0;

%--
% thermal velocity
vt = sqrt((T_ev*e) / me);

%--
% ion constants (95% D, 5% H)
qd = abs(e);
mp = 1.67e-27;
md = 2.0*mp;

qh = abs(e);
mh = mp;

%-- 
% cyclotron frequencies
om_ce = e*B0/me;
om_cd = qd*B0/md;
om_ch = qh*B0/mh;

%--
% rotation matrix
alpha = 0.5;
beta = 0.5;

r11 = cos(beta)*cos(alpha);
r12 = cos(beta)*sin(alpha);
r13 = -sin(beta);
r21 = -sin(alpha);
r22 = cos(alpha);
r23 = 0.0;
r31 = sin(beta)*cos(alpha);
r32 = sin(beta)*sin(alpha);
r33 = cos(beta);

r = [[r11, r12, r13]
     [r21, r22, r23]
     [r31, r32, r33]];
 
r_para = r(3,:);
Bvec = B0*r_para;

%--
% electron calcs; density, plasma frequency
Ne = N0;
om_pe = sqrt(Ne*e^2/(me*eps0));

%--
% ion calcs (95% D, 5% H); density, plasma frequency
Nd = 0.95*N0;
om_pd = sqrt(Nd*qd^2/(md*eps0));
Nh = 0.05*N0;
om_ph = sqrt(Nh*qh^2/(mh*eps0));

%--
% poisson solve for static potential -- solution "static_pot"

poisson_sol

%--
% static electric field calculation



%--
% wave solver to find rf electric field

wave_sol