%----------------------------------------%
% main script for 1D DVE loop            %
% rlbarnett c3149416 061117              %
%----------------------------------------%

%clear all; close all

%--
% constants
mu0 = 4.0*pi*1.0e-7;
eps0 = 8.85e-12;
c0 = 1.0/sqrt(eps0*mu0);

%--
% magnetic field (tesla)
B0 = 2.6;

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
lamby = 1.0;
lambz = 1.0;

%--
% spatial domain
npts = 128;
dx = 0.002;
xmin = 0.0;
xmax = 0.2;
% npts = ((xmax - xmin)/dx);
xax = linspace(xmin, xmax, npts);
% xax = xmin:dx:xmax;
% dx = (xmax - xmin)/(npts - 1);

%--
% background density -- set to zero for vacuum case
Nmax = 5.0e17;
N0 = Nmax*ones(npts,1);

%--
% initialise perturbed density as zero
N1 = zeros(npts,1);
N1e = N1;
N1i = N1;

%--
% initialise perturbed velocity as zero
v1 = zeros(npts,1);
v1e = v1;
v1i = v1;

%--
% electron constants
e = -1.6022e-19;
me = 9.11e-31;

%--
% temperature
T_ev = 15.0;

%--
% thermal velocity
vt = sqrt((T_ev*abs(e)) / me);

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
% temporal step
period = 2.0*pi / om_ce;
t = 0.0;
num_cyc = 10;
num_points = 100;
dt = period/num_points;
tmax = num_cyc*num_points*dt;
nmax = ((tmax - t) / dt);

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

rot = [[r11, r12, r13]
     [r21, r22, r23]
     [r31, r32, r33]];
 
e_para = rot(3,:);
Bvec = B0*e_para;

%--
% electron calcs; density, plasma frequency 
N0e = 1.05*N0;
om_pe = sqrt(N0e*e^2/(me*eps0));

%--
% ion calcs (95% D, 5% H); density, plasma frequency
N0d = 0.95*N0;
om_pd = sqrt(N0d*qd^2/(md*eps0));
N0h = 0.05*N0;
om_ph = sqrt(N0h*qh^2/(mh*eps0));
N0i = N0h + N0d;

%--
% poisson solve for static potential -- solution (output) "static_pot"

poisson_sol;%(Ne, Nh, Nd, lamby, lambz, e, eps0, dx, npts);

%--
% static electric field calculation -- solution (output) "static_e(x,y,z)"

static_e;

%--
% wave solver to find rf electric field -- solution (output) "rf_e(x,y,z)"

wave_sol;

%--
% ponderomotive acceleration calculation -- solution (output)
% "a_pond(x,y,z)"

a_pond;

%--
% pressure term -- solution (output) "press(x,y,z)"

pressure;

%--
% calculate perpendicular drift velocities analytically -- solution
% (output) "v_perp(1,2)"

v_drift;

%--
% solve equation 23 (DVE 2015): slow time scale continuity equation yielding v
% parallel -- solution (output)

cont_slow;