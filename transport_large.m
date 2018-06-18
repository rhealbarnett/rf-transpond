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
nu = 1000.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 1.0;

%-- include two additional gridpoints for the density ghost points
%-- velocity grid will then be defined as having npts-1 (xax(1:npts-1))
%-- density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 128;
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
Nmax = 17;
% Nmin = 16;
% slope = (Nmax - Nmin) ./ (xmax - xmin);
% n_new = (10^Nmax - 10^Nmin)*exp(-10.0*nxax) + 10^Nmin;
n_new = (10^Nmax)*ones(1,npts);
dnx = gradient(n_new,nxax);

%-- density source
rate_coeff = 10e-14;
rate_min = 10^0.0;
rate_max = (10^Nmax);
% n_neut = (rate_max - rate_min)*exp(-90.0*nxax(1,1:end/2)) + rate_min;
n_neut = zeros(1,npts);
n_neut(1:round(2*npts/5)+1) = 10.^(-42.5*nxax(1:round(2*npts/5)+1) + Nmax);
n_neut(end-round(2*npts/5):end) = fliplr(n_neut(1:round(2*npts/5)+1));
% n_neut = [n_neut,fliplr(n_neut)];
n_neut = n_neut';
n_source = zeros(npts,1);

for ii=1:npts
    n_source(ii,1) = n_neut(ii,1)*n_new(1,ii)*rate_coeff;
end

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
tmax = 1.0e-2;
if (0.99*(dx^2)/(2.0*nu))<(0.99*dx/max(abs(vx_new)))
    dt = 0.99*(dx^2)/(2.0*nu);
elseif (0.99*(dx^2)/(2.0*nu))>(0.99*dx/max(abs(vx_new)))
    dt = 0.99*dx/max(abs(vx_new));
end
nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);

