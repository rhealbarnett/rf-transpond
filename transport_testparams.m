%-------------------------------------------------------%
% parameter file for coupled transport equations        %
% test case:                                            %
% domain inversely proportional to cs                   %
% constant velocity -- zero                             %
% density 'blob' -- 1 <= n <= 2                         %
% rlbarnett c3149416 140518                             %
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
cs = 10*sqrt((Te + Ti)*e/m);
nu = 1.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 10/cs;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% total 128 grid points
% density solution space will be defined as having npts-2 (xax(2:npts-1))
% -- total 127 gridpoints
npts = 1025;
dx = (xmax - xmin)/(npts - 1);
% xax = linspace(xmin,xmax,npts);
nxax = linspace(xmin-0.5*dx,xmax+0.5*dx,npts);
vxax = linspace(xmin,xmax,npts-1);

%------
% temporal domain %
%------
tmin = 0;

% set dt based on CFL conditions, check during loop if violated
nmax = 1.0e5;
dt = 0.99*(dx^2)/(2.0*nu);

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
%-------------------------------------------------------------------------%

% initial density profile -- Gaussian 'blob' at centre of the domain
n_new = normpdf(nxax,(xmax+0.5*dx)/2,(xmax+0.5*dx)/30);
n_new = n_new/max(n_new);

% shift -- density used in denominator for calculation, can't have zeros 
n_new = n_new + 1.0;
dnx = gradient(n_new,nxax);

% initial velocity
vx_new = zeros(1,npts-1);
% vx_ax = linspace(0,1,npts-1);
% vx_new = -(cs/2)*vx_ax + cs;
% vx_new = cs*ones(1,npts-1);

% initialise coefficient matrices for density, velocity, and momentum equation 
% rhs 'source' term
nA = zeros(npts,npts);
vxA = zeros(npts-1,npts-1);
vx_source = zeros(npts-1,1);

% fill boundary conditions in coefficient matrix
% zero flux on density ghost points
nA(1,1) = 1.0;
nA(1,2) = -1.0;
nA(end,end) = 1.0;
nA(end,end-1) = -1.0;
% Dirichlet conditions on velocity 
vxA(1,1) = 1.0;
vxA(end,end) = 1.0;

