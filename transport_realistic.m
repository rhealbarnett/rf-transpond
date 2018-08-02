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
cs = sqrt((Te + Ti)*e/m);
nu = 0.0;%1000.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 0.1;

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
% N_grad = (Nmax - Nmin)/(xmax - xmin);
% slope = (Nmax - Nmin) ./ (xmax - xmin);
% n_new = (N_grad)*exp(10*nxax(2:npts-1));
% n_new = 10.^(N_grad*nxax(2:npts-1) + Nmax);
n_new = (Nmax)*ones(1,npts-2);
dnx = gradient(n_new,nxax(2:npts-1));

%-- density source
rate_coeff = 10e-14;
decay_index = round((npts-2)/2.5);
cosax = linspace(0,pi,decay_index+1);
neut_max = 17.5;
neut_min = 14;
decay_length = 0.4;
decay_gradient = (neut_min - neut_max)/decay_length;
% n_neut = (rate_max - rate_min)*exp(-90.0*nxax(1,1:end/2)) + rate_min;
n_neut = zeros(1,npts-2);
% n_neut(1:decay_index + 1) = 10.^(decay_gradient*nxax(1:decay_index + 1) + neut_max);
n_neut(1:decay_index+1) = 10^neut_max*(cos(cosax)+1.01)/2;%.*exp(-4*cosax);
n_neut(end-decay_index:end) = fliplr(n_neut(1:decay_index + 1));
% n_neut = [n_neut,fliplr(n_neut)];
n_neut = n_neut';
n_neut(1:end-decay_index) = n_neut(decay_index)/2;
n_source = zeros((npts-2),1);

% for ii=1:npts-2
%     n_source(ii,1) = n_neut(ii,1)*n_neut(ii,1)*rate_coeff;
% end

%-- initial velocity
% vx_ax = linspace(-1,0,npts-1);
% vx_new = (cs)*vx_ax;
% vx_new = cs*ones(1,npts-1);
vx_new = 400*cs*vxax.^2 - 40*cs*vxax + cs;
% vx_new = 400*cs*vxax.^2 - 40*cs*vxax;
% vx_new = -400*cs*vxax.^2 + 40*cs*vxax - cs;
vx_new(1,1) = cs;
% vx_new(1,end) = cs;

%-- initialise coefficient matrices for density, velocity, and momentum equation 
%-- rhs 'source' term
nA = zeros(npts,npts);
vxA = zeros(npts-1,npts-1);
vx_source = zeros(npts-1,1);

%-- fill boundary conditions in coefficient matrix
%-- Dirichlet conditions on velocity 
vxA(1,1) = 1.0;
% vxA(end,end) = 1.0;

%-- set dt based on CFL conditions, check during loop if violated
tmax = 1.0e-6;
if (0.8*(dx^2)/(2.0*nu))<(0.8*dx/max(abs(vx_new)))
    dt = 0.8*(dx^2)/(2.0*nu);
elseif (0.8*(dx^2)/(2.0*nu))>(0.8*dx/max(abs(vx_new)))
    dt = 0.8*dx/max(abs(vx_new));
end
nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);

%%

%-----------------------------------------------------------------------------%
% Set up electric field profile & bits for PF                                 %
% Using E||,max = 300Vcm^-1 from J. Myra 2006 Nonlinear interactions paper    %
% exponential decay away from antenna location                                %
%-----------------------------------------------------------------------------%

Emax = 3.0e5;
freq = 50.0e6;
om = 2.0*pi*freq;

% Efield = exp(1.0e3*vxax);
Efield = (Emax/2)*(-cos(cosax)+1.01);
Efield = [(zeros(1,npts-1-length(cosax))), Efield];
Efield = Efield./max(Efield);
Efield = Emax*Efield;
Efield = Efield.^2;
Efield = zeros(npts-1);

pond_const = (1.0/4.0)*((e^2)/(m*om^2));


