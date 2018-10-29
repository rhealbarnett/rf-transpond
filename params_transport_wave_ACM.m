%-----------------------------------------------------------------------%
% parameter file for coupled transport equations & wave eq        
% 'realistic' case:                                     
% domain ~cms                                           
% velocity and density taken from equilibrium solve of transport eqs
% parameters for Alcator Cmod
% rlbarnett c3149416 291018                             
%-----------------------------------------------------------------------%

%------
% constants %
%------
constants;
m = const.mp;
e = const.e;
c0 = const.c0;
eps0 = const.eps0;

%------
% transport parameters %
%------
Te = 5.0;
Ti = 5.0;
cs = sqrt((Te + Ti)*e/m);
nu = 1.0;

%------
% wave solver parameters %
%------
freq = 4.6e9;
om = 2*pi*freq;
k0 = om/c0;
k_para = 183;
B0 = 4.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 0.1;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 4096;
dx = (xmax - xmin)/(npts - 1);
nxax = linspace(xmin-0.5*dx,xmax+0.5*dx,npts);
vxax = linspace(xmin,xmax,npts-1);

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
%-------------------------------------------------------------------------%

plot_num = 10;
equib = load('equib.mat');

%-- initial density profile
n_new = equib.n(plot_num+2,:);
n_new = full(n_new)*10;
n_init = n_new;

%-- initial velocity
vx_new = equib.vx(plot_num+2,:);
vx_new = full(vx_new);
vx_init = vx_new;

%-- density source
rate_coeff = 1.*10^-14;
decay_index = round(npts/4);
cosax = linspace(pi,2*pi,decay_index);
neut_max = 0.5e18;
neut_min = 1.0e14;
decay_length = 0.4;
decay_gradient = (neut_min - neut_max)/decay_length;
n_neut = zeros(1,npts);
n_neut(end-decay_index+1:end) = neut_max*(cos(cosax) + 1)/2;
n_neut(1:end-decay_index+1) = n_neut(end-decay_index+2);
n_source = zeros(1,npts);

for jj=2:npts-1
    n_source(jj) = n_new(jj)*n_neut(jj)*rate_coeff;
end

ns_mult = n_neut(end-1)/n_new(end-1);
n_source = (n_source*ns_mult*150);

%%
%-- initialise coefficient matrices for density, velocity, and momentum equation 
%-- rhs 'source' term
nA = sparse(npts,npts);
nI = sparse(eye(npts,npts));
vx_pos = sparse(npts-1,npts-1);
vx_neg = sparse(npts-1,npts-1);
vx_diff = sparse(npts-1,npts-1);
vx_I = sparse(eye(npts-1,npts-1));

%-- set dt based on CFL conditions, check during loop if violated
tmin = 0.0;
tmax = 4.0e-6;
cfl_fact = 0.99;

if ((cfl_fact*(dx^2)/(2.0*nu))<(cfl_fact*dx/max(abs(vx_new))))
    dt = cfl_fact*(dx^2)/(2.0*nu);
elseif (cfl_fact*(dx^2)/(2.0*nu))>(cfl_fact*dx/max(abs(vx_new)))
    dt = cfl_fact*dx/max(abs(vx_new));
else
    dt = cfl_fact*dx/cs;
end

dt = 2.0*dt;
 
nmax = round(tmax/dt);
plot_freq = int64(nmax/plot_num);
tax = linspace(tmin,tmax,nmax);
mult = 1.0/dx;

%%

%-----------------------------------------------------------------------------%
% Set up electric field profile & bits for PF                                 %
% Using E||,max = 300Vcm^-1 from J. Myra 2006 Nonlinear interactions paper    %
% exponential decay away from antenna location                                %
%-----------------------------------------------------------------------------%

Emax = 3.0e4;

Efield = exp(100*vxax);
Efield = Efield./max(Efield);
Efield = Emax*Efield;
Efield = Efield.^2;
% ----------- %
% set e field to zero for testing 
% ----------- %
% Efield = zeros(1,npts-1);

pond_const = (1.0/4.0)*((e^2)/(m*om^2));

vx_mat = sparse(nmax,npts-1);
n_mat = sparse(nmax,npts);
pressure_mat = sparse(nmax,npts-2);

vx_mat(1,:) = vx_new;
n_mat(1,:) = n_new;



