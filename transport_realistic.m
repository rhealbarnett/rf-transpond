%-------------------------------------------------------%
% parameter file for coupled transport equations        %
% 'realistic' case:                                     %
% domain ~cms                                           %
% linear velocity -- cs/2 <= v <= cs                    %
% exponential density profile -- 1e16 <= n <= 1e17      %
% rlbarnett c3149416 150518                             %
%-------------------------------------------------------%

%------
% constants %
%------
constants;
m = mp;

%------
% parameters %
%------
Te = 10.0;
Ti = 5.0;
T = Te + Ti;
% cs = sqrt((Te + Ti)*e/m);
cs = 10;
% nu = 1000.0;
nu = 0.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 1.0;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 32;
dx = (xmax - xmin)/(npts - 1);
nxax = linspace(xmin-0.5*dx,xmax+0.5*dx,npts);
vxax = linspace(xmin,xmax,npts-1);

%------
% temporal domain %
%------
tmin = 0;

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
%-------------------------------------------------------------------------%

%-- initial density profile
Nmax = 16;
Nmin = 15;
n_new = (10^Nmax)*ones(1,npts);


%-- density source
rate_coeff = 10e-14;
decay_index = round((npts));
cosax = linspace(0,2*pi,decay_index);
neut_max = 13.2;
neut_min = 14;
decay_length = 0.4;
decay_gradient = (neut_min - neut_max)/decay_length;
n_neut = zeros(1,npts);
n_neut(1:decay_index) = 10^neut_max*(cos(cosax)+1.01)/2;%.*exp(-4*cosax);
% n_neut(end-decay_index:end) = fliplr(n_neut(1:decay_index + 1));
% n_neut(decay_index+1:end-decay_index) = (n_neut(decay_index)/2);
n_source = zeros((npts),1);

for ii=1:npts
    n_source(ii,1) = n_neut(1,ii)*n_neut(1,ii)*rate_coeff;
end

%-- initial velocity
vx_ax = linspace(-1,1,npts-1);
vx_new = (cs)*vx_ax;

%-- initialise coefficient matrices for density, velocity, and momentum equation 
%-- rhs 'source' term
nA = zeros(npts,npts);
vxA = zeros(npts-1,npts-1);
vx_source = zeros(npts-1,1);


%-- set dt based on CFL conditions, check during loop if violated
tmax = 2.0e-4;
cfl_fact = 0.99;

if ((cfl_fact*(dx^2)/(2.0*nu))<(cfl_fact*dx/max(abs(vx_new))))
    dt = cfl_fact*(dx^2)/(2.0*nu);
elseif (cfl_fact*(dx^2)/(2.0*nu))>(cfl_fact*dx/max(abs(vx_new)))
    dt = cfl_fact*dx/max(abs(vx_new));
else
    dt = cfl_fact*dx/cs;
end

nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);
mult = dt/dx;

%%

%-----------------------------------------------------------------------------%
% Set up electric field profile & bits for PF                                 %
% Using E||,max = 300Vcm^-1 from J. Myra 2006 Nonlinear interactions paper    %
% exponential decay away from antenna location                                %
%-----------------------------------------------------------------------------%

Emax = 3.0e4;
freq = 50.0e6;
om = 2.0*pi*freq;

% Efield = exp(1.0e3*vxax);
Efield = (Emax/2)*(-cos(cosax)+1.01);
Efield = [(zeros(1,npts-1-length(cosax))), Efield];
Efield = Efield./max(Efield);
Efield = Emax*Efield;
Efield = Efield.^2;
% ----------- %
% set e field to zero for testing 
% ----------- %
Efield = zeros(1,npts-1);

% pond_const = (1.0/4.0)*((e^2)/(m*om^2));
% pond_source = zeros(npts-1,1);


