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
nu = 10.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 1.0;

%-- include two additional gridpoints for the density ghost points
%-- velocity grid will then be defined as having npts-1 (xax(1:npts-1))
%-- density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 4072;
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

%-- initial density profile -- Gaussian 'blob' at centre of the domain
Nmax = 17;
% Nmin = 16;
% slope = (Nmax - Nmin) ./ (xmax - xmin);
% n_new = (10^Nmax - 10^Nmin)*exp(-10.0*nxax) + 10^Nmin;
n_new = (10^Nmax)*ones(1,npts);
dnx = gradient(n_new,nxax);

%-- density source
rate_coeff = 10e-14;
rate_min = 10^0.0;
rate_max = (10^Nmax)^2*rate_coeff;
n_source = (rate_max - rate_min)*exp(-90.0*nxax(1,1:end/2)) + rate_min;
% n_source = rate_max*(1.0e-2).^nxax(1:end/2);
n_source = [n_source,fliplr(n_source)];
n_source = n_source';
% n_source = zeros(npts,1);


%-- initial velocity
% vx_ax = linspace(0,1,npts-1);
% vx_new = (2.0*cs)*vx_ax - cs;
% vx_new = cs/2*ones(1,npts-1);
vx_new = zeros(1,npts-1);
vx_new(1,1) = -cs;
vx_new(1,end) = cs;

%-- initialise coefficient matrices for density, velocity, and momentum equation 
%-- rhs 'source' term
nA = zeros(npts,npts);
vxA = zeros(npts-1,npts-1);
vx_source = zeros(npts-1,1);

%-- fill boundary conditions in coefficient matrix
%-- zero flux on density ghost point at xmax
%-- Dirichlet at xmin
nA(1,1) = 1.0;
nA(1,2) = -1.0;
nA(end,end) = 1.0;
nA(end,end-1) = -1.0;
%-- Dirichlet conditions on velocity 
vxA(1,1) = 1.0;
vxA(end,end) = 1.0;

%-- set dt based on CFL conditions, check during loop if violated
tmax = 1.0e-4;
if (0.99*(dx^2)/(2.0*nu))<(0.99*dx/max(abs(vx_new)))
    dt = 0.99*(dx^2)/(2.0*nu);
elseif (0.99*(dx^2)/(2.0*nu))>(0.99*dx/max(abs(vx_new)))
    dt = 0.99*dx/max(abs(vx_new));
end
nmax = round(tmax/dt);

