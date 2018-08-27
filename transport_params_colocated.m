%-------------------------------------------------------%
% parameter file for coupled transport equations        %
%                                                       %
% domain ~cms                                           %
% linear velocity -- cs/2 <= v <= cs                    %
% flat density profile                                  %
% rlbarnett c3149416 130818                             %
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
% cs = 10;
nu = 1000.0;
% nu = 0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 1.0;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 128;
dx = (xmax - xmin)/(npts - 1);
nxax = linspace(xmin,xmax,npts-1);
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
% n_new = normpdf(nxax,(xmax)/2,(xmax)/30);
% n_new = (10^Nmax)*n_new/max(n_new);
n_new = (10^Nmax)*ones(1,npts-1);

%-- density source
rate_coeff = 10e-14;
decay_index = round((npts-1)/2.5);
cosax = linspace(0,pi,decay_index+1);
neut_max = 16.0;
neut_min = 14;
decay_length = 0.4;
decay_gradient = (neut_min - neut_max)/decay_length;
% n_neut = (rate_max - rate_min)*exp(-90.0*nxax(1,1:end/2)) + rate_min;
n_neut = zeros(1,npts-1);
% n_neut(1:decay_index + 1) = 10.^(decay_gradient*nxax(1:decay_index + 1) + neut_max);
n_neut(1:decay_index+1) = (10^neut_max)*(cos(cosax)+1.01)/2;
n_neut(end-decay_index:end) = fliplr(n_neut(1:decay_index + 1));
% % n_neut = [n_neut,fliplr(n_neut)];
% n_neut = n_neut';
n_neut(decay_index+2:end-decay_index-1) = (n_neut(decay_index)/2);
% n_neut = fliplr(n_neut);
n_source = zeros((npts-1),1);

% for ii=1:npts-1
%     n_source(ii,1) = n_neut(1,ii)*n_neut(1,ii)*rate_coeff;
% end

%-- initial velocity
vx_ax = linspace(-1,1,npts-1);
vx_new = cs*vx_ax;
% vx_new(1,1) = -cs/2;
% vx_new = zeros(1,npts-1);
% vx_new(1:(npts-1)/2) = -cs;
% vx_new((npts-1)/2:npts-1) = cs;


%-- initialise coefficient matrices for density, velocity, and momentum equation 
%-- rhs 'source' term
nA = zeros(npts-1,npts-1);
vxA = zeros(npts-1,npts-1);
vx_source = zeros(npts-1,1);

%-- set dt based on CFL conditions, check during loop if violated
tmax = 1.0e-5;
cfl_fact = 0.8;

if ((cfl_fact*(dx^2)/(2.0*nu))<(cfl_fact*dx/max(abs(vx_new))))
    dt = cfl_fact*(dx^2)/(2.0*nu);
elseif (cfl_fact*(dx^2)/(2.0*nu))>(cfl_fact*dx/max(abs(vx_new)))
    dt = cfl_fact*dx/max(abs(vx_new));
else
    dt = cfl_fact*dx/cs;
end

mult = dt/dx;
nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);

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
Efield = zeros(1,npts);

pond_const = (1.0/4.0)*((e^2)/(m*om^2));
pond_source = zeros(npts,1);


