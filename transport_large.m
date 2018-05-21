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
nu = 100.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 1.0;

%-- include two additional gridpoints for the density ghost points
%-- velocity grid will then be defined as having npts-1 (xax(1:npts-1))
%-- density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 8142;
dx = (xmax - xmin)/(npts - 1);
% r = 0.99;
% dx_arr = zeros(1,npts);
% dx_arr(1,1) = dx;
% for ii=2:npts
%     dx_arr(1,ii) = r*dx_arr(1,ii-1);
% end
% % xax = linspace(xmin,xmax,npts);
% nxax = zeros(1,npts);
% dx_min = dx_arr(1,1) - r*dx_arr(1,1);
% dx_max = dx_arr(1,npts) + r*dx_arr(1,npts);
% nxax(1,1) = xmin - 0.5*dx_min;
% nxax(1,npts) = xmax + 0.5*dx_max;
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
tmax = 1.0e-6;
if (0.99*(dx^2)/(2.0*nu))<(0.99*dx/max(abs(vx_new)))
    dt = 0.99*(dx^2)/(2.0*nu);
elseif (0.99*(dx^2)/(2.0*nu))>(0.99*dx/max(abs(vx_new)))
    dt = 0.99*dx/max(abs(vx_new));
end
nmax = round(tmax/dt);

