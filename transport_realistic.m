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
% parameters
%------
Te = 10.0;
Ti = 5.0;
T = Te + Ti;
cs = sqrt((Te + Ti)*e/m);
nu = 1.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 0.01;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 512;
dx = (xmax - xmin)/(npts - 1);
% xax = linspace(xmin,xmax,npts);
nxax = linspace(xmin-0.5*dx,xmax+0.5*dx,npts);
vxax = linspace(xmin,xmax,npts-1);

%------
% temporal domain %
%------
tmin = 0;

% set dt based on CFL conditions, check during loop if violated
tmax = 1.0e-3;
dt = 0.99*(dx^2)/(2.0*nu);
nmax = int64(tmax/dt);
% nmax = 1.0e5;

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
%-------------------------------------------------------------------------%

% initial density profile -- Gaussian 'blob' at centre of the domain
Nmax = 17;
Nmin = 16;
slope = (Nmax - Nmin) ./ (xmax - xmin);
n_new = 10.^(-slope*nxax + Nmax);
dnx = gradient(n_new,nxax);

% initial velocity
% vx_ax = linspace(0,1,npts-1);
% vx_new = (cs/2)*vx_ax + cs/2;
vx_new = cs/2*ones(1,npts-1);
vx_new(1,end) = cs;

% initialise coefficient matrices for density, velocity, and momentum equation 
% rhs 'source' term
nA = zeros(npts,npts);
vxA = zeros(npts-1,npts-1);
vx_source = zeros(npts-1,1);

% fill boundary conditions in coefficient matrix
% zero flux on density ghost point at xmax
% Dirichlet at xmin
nA(1,1) = 1.0;
nA(end,end) = 1.0;
nA(end,end-1) = -1.0;
% Dirichlet conditions on velocity 
vxA(1,1) = 1.0;
vxA(end,end) = 1.0;

