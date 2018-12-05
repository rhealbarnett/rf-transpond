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


% function [n_init, vx_init] = transport_test(params,xmin,xmax,npts)

const = constants();
m = const.mp;
e = const.e;

%------
% parameters %
%------
% params = transport_params();
% Te = params.Te;
% Ti = params.Ti;
% nu = params.nu;
% cs = params.cs;

Te = 5.0;
Ti = 10.0;
nu = 1.0;
cs = sqrt((Te + Ti)*e/m);

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
% nxax = linspace(xmin,xmax,npts-1);
vxax = linspace(xmin,xmax,npts-1);

%%

%------
% temporal domain %
%------
tmin = 0;

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
%-------------------------------------------------------------------------%

% equib = load('equib2.mat');

%-- initial density profile
Nmax = (1.0e18);
Nmin = (0.5e18);
% n_new = (Nmin*(fliplr(nxax)/max(nxax)) + Nmin);
n_new = Nmax*(1.0 - (1.0e-5)*exp((nxax/max(nxax))*10));
% n_new = (Nmax*ones(1,npts));
% n_new = equib.n_new;
% n_new = interp1(nxax(2:npts-1),n_new(2:npts-1),nxax,'linear','extrap');
% n_new = full(n_new)*10;
n_avg = interp1(nxax,n_new,vxax);
n_init = n_new;

%-- initial velocity
% vx_new = (cs)*(vxax/max(vxax));
vx_new = 1.0 + exp(vxax*100);
vx_new = cs*(vx_new/max(vx_new));
% vx_new = zeros(1,npts-1);
% vx_new(1,end) = cs;
% vx_new = equib.vx_new;
% vx_new = full(vx_new);
vx_init = vx_new;


%%
%-------------------------------------------------------------------------%
% CALCULATE DENSITY SOURCE                                                %
%-------------------------------------------------------------------------%

% rate coefficient (constant)
rate_coeff = (1.0e-14);
% approx size of non-zero portion of neutral profile (1/4 domain)
decay_index = round((npts)/4);
% calculate shape of neutral profile
cosax = linspace(pi,2*pi,decay_index);
% max neutral value (at wall)
neut_max = (1.0e18);

% initialise and fill neutral density array
n_neut = zeros(1,npts);
n_neut(end-decay_index+1:end) = neut_max*((cos(cosax) + 1)/2);
n_neut(1:end-decay_index+1) = n_neut(end-decay_index+2);

% calculate density source
n_source = (n_new.*(n_neut)*(rate_coeff));

% interpolate source onto velocity grid
source_avg = interp1(nxax,n_source,vxax);
% calculate integral of density source over the velocity grid
source_int = trapz(vxax,source_avg);
% source_int = trapz(nxax,n_source);
% calculate flux at the rh boundary (wall)
% rflux = n_avg(end)*vx_new(end);
rflux = n_avg(end)*vx_new(end);
% calculate the constant multiplier to match density out = in
ns_mult = rflux/source_int;
% multiply n0(x)n(x,t) by the constant calculated in previous step
n_source = (n_source*ns_mult)*0.5;
nv_source = source_avg*ns_mult;
n_source(1,1) = 0.0; n_source(1,end) = 0.0;

% n_source = equib.n_source;

fprintf('Initial outward flux at RH boundary %d\n',rflux)
fprintf('Initial density source integral %d\n',source_int)
fprintf('Density source integral after normalisation %e\n',trapz(nxax,n_source))
fprintf('Initial total number of particles %e\n',trapz(nxax,n_new))
fprintf(['Difference between outward flux at RH boundary and integral ' ...
    ' of density source %e\n'], rflux - trapz(vxax,nv_source))    


%%
%-------------------------------------------------------------------------%
% Initialise coefficient matrices                                         %
%-------------------------------------------------------------------------%

%-- initialise coefficient matrices as sparse
nA = sparse(npts,npts);
nI = sparse(eye(npts,npts));
vx_pos = sparse(npts-1,npts-1);
vx_neg = sparse(npts-1,npts-1);
vx_diff = sparse(npts-1,npts-1);
vx_I = sparse(eye(npts-1,npts-1));

%-------------------------------------------------------------------------%
% Calculate time step                                                                  %
%-------------------------------------------------------------------------%
%-- set dt based on CFL conditions, check during loop if violated
tmax = 2.0e-5;
cfl_fact = 0.99;

if ((cfl_fact*((dx)^2)/(2.0*nu))<(cfl_fact*(dx)/max(abs(vx_new))))
    dt = cfl_fact*((dx)^2)/(2.0*nu);
elseif (cfl_fact*(dx)^2/(2.0*nu))>(cfl_fact*(dx)/max(abs(vx_new)))
    dt = cfl_fact*(dx)/max(abs(vx_new));
else
    dt = cfl_fact*(dx)/cs;
end

dt = 2.0*dt;
nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);
mult = 1.0/dx;

%%

%-----------------------------------------------------------------------------%
% Set up electric field profile & bits for PF                                 %
% Using E||,max = 300Vcm^-1 from J. Myra 2006 Nonlinear interactions paper    %
% exponential decay away from antenna location                                %
%-----------------------------------------------------------------------------%

Emax = 6.0e4;
freq = 80.0e6;
om = 2.0*pi*freq;

% Efield = exp(1.0e3*vxax);
% Efield = (Emax/2)*(cos(cosax)+1.01);
% Efield = [(zeros(1,npts-1-length(cosax))), Efield];
% Efield = exp(100*vxax);
% Efield = Efield./max(Efield);
% Efield = Emax*Efield;
% Efield = interp1(nxax,real(rf_ex),vxax);
% Efield = Efield.^2;
% ----------- %
% set e field to zero for testing 
% ----------- %
Efield = zeros(1,npts-1);

pond_const = (1.0/4.0)*((e^2)/(m*om^2));
% pond_source = zeros(1,npts-1);

vx_mat = sparse(nmax,npts-1);
n_mat = sparse(nmax,npts);
pressure_mat = sparse(nmax,npts-2);

vx_mat(1,:) = vx_new;
n_mat(1,:) = n_new;

% end

