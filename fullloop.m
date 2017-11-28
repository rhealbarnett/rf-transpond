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
npts = 1000;
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

%--
% solve equation 24 (DVE 2015): slow time scale parallel equation of motion yielding log(N0) 
% -- solution (output)

% eqofmot_slow;

%--
% figure to follow important quantities as changes are made
figure(1);
set(gcf,'Position',get(0,'Screensize'))

subplot(2,3,1)
plot(xax,static_pot,'r')
ylabel('$\phi$ (V)','Fontsize',16)
ytickformat('%.2f')

subplot(2,3,2)
plot(xax,N0e,'k')
ylabel('N$_{0,e}$','Fontsize',16)
ytickformat('%.2f')

subplot(2,3,3)
plot(xax,N1e,'k')
ylabel('N$_{1,e}$','Fontsize',16)
ytickformat('%.2f')

subplot(2,3,4)
plot(xax,static_ex,'k')
ylabel('E$_{0,ex}$','Fontsize',16)
ytickformat('%.2f')

subplot(2,3,5)
plot(xax,static_ey,'k')
ylabel('E$_{0,ey}$','Fontsize',16)
ytickformat('%.2f')

subplot(2,3,6)
plot(xax,static_ez,'k')
ylabel('E$_{0,ez}$','Fontsize',16)
ytickformat('%.2f')

figure(2);
set(gcf,'Position',get(0,'Screensize'))

subplot(3,3,1)
plot(xax,real(rf_ex),'k')
hold on
plot(xax,imag(rf_ex),'--b')
legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest') 
ylabel('E$_{1,ex}$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,2)
plot(xax,real(rf_ey),'k')
hold on
plot(xax,imag(rf_ey),'--b')
ylabel('E$_{1,ey}$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,3)
plot(xax,real(rf_ez),'k')
hold on
plot(xax,imag(rf_ez),'--b')
ylabel('E$_{1,ez}$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,4)
plot(xax,pond_pote,'r')
ylabel('$\Theta$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,7)
plot(xax,a_pondex,'k')
ylabel('$(-\nabla\Theta)_x$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,8)
plot(xax,a_pondey,'k')
ylabel('$(-\nabla\Theta)_y$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,9)
plot(xax,a_pondez,'k')
ylabel('$(-\nabla\Theta)_z$','Fontsize',16)
ytickformat('%.2f')

figure(3);
set(gcf,'Position',get(0,'Screensize'))

subplot(3,3,1)
plot(xax,vd_perp1e,'k')
ylabel('v$_{e\perp,1}$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,2)
plot(xax,vd_perp2e,'k')
ylabel('v$_{e\perp,2e}$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,3)
plot(xax,v_parae,'k')
ylabel('v$_{e\parallel}$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,4)
plot(xax,v1e,'r')
ylabel('v$_{1,e}$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,7)
plot(xax,pressex,'k')
ylabel('$(\nabla N_0/N_0)_x$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,8)
plot(xax,pressey,'k')
ylabel('$(\nabla N_0/N_0)_y$','Fontsize',16)
ytickformat('%.2f')

subplot(3,3,9)
plot(xax,pressez,'k')
ylabel('$(\nabla N_0/N_0)_z$','Fontsize',16)
ytickformat('%.2f')




