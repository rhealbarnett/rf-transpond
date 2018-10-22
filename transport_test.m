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
m = const.mp;
e = const.e;

%------
% parameters %
%------
Te = 10.0;
Ti = 5.0;
cs = sqrt((Te + Ti)*e/m);
nu = 1.0;

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
Nmax = 1.0e18;
Nmin = 0.5e18;
n_new = Nmin*(fliplr(nxax)/max(nxax)) + Nmin;
n_avg = (n_new(1:npts-1) + n_new(2:npts))/2.0;

%-- initial velocity
vx_new = (cs)*(vxax/max(vxax));
% vx_new = zeros(1,npts-1);
% vx_new(1,end) = cs;

%-- flux at boundaries
fl = vx_new(1,1)*((n_new(1,1)+n_new(1,2))/2);
fr = vx_new(1,end)*((n_new(1,end) + n_new(1,end-1))/2);
ft = fr - fl;
flux = vx_new.*n_avg;

%-- density source
rate_coeff = 10^-14;
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

source_int = trapz(n_source);

flux_int = trapz(flux);
ns_mult = n_neut(end-1)/n_new(end-1);
n_source = n_source*ns_mult*15;

% if source_int~=ft
%     diff = ft - source_int;
%     bal = diff/(npts-2);
%     source_bal = bal*ones(1,npts-2);
%     n_source(2:npts-1) = n_source(2:npts-1) + source_bal;
% end
    

%-- initialise coefficient matrices for density, velocity, and momentum equation 
%-- rhs 'source' term
nA = sparse(npts,npts);
nI = sparse(eye(npts,npts));
vx_pos = sparse(npts-1,npts-1);
vx_neg = sparse(npts-1,npts-1);
vx_diff = sparse(npts-1,npts-1);
vx_I = sparse(eye(npts-1,npts-1));

%-- set dt based on CFL conditions, check during loop if violated
tmax = 1.0e-4;
cfl_fact = 0.99;

if ((cfl_fact*(dx^2)/(2.0*nu))<(cfl_fact*dx/max(abs(vx_new))))
    dt = cfl_fact*(dx^2)/(2.0*nu);
elseif (cfl_fact*(dx^2)/(2.0*nu))>(cfl_fact*dx/max(abs(vx_new)))
    dt = cfl_fact*dx/max(abs(vx_new));
else
    dt = cfl_fact*dx/cs;
end

dt = 2.0*dt;
% nmax = 20400;   
nmax = round(tmax/dt);
% tmax = nmax*dt;
tax = linspace(tmin,tmax,nmax);
mult = 1.0/dx;
tol = 1.0e-2;

%%

%-----------------------------------------------------------------------------%
% Set up electric field profile & bits for PF                                 %
% Using E||,max = 300Vcm^-1 from J. Myra 2006 Nonlinear interactions paper    %
% exponential decay away from antenna location                                %
%-----------------------------------------------------------------------------%

Emax = 6.0e4;
freq = 50.0e6;
om = 2.0*pi*freq;

% Efield = exp(1.0e3*vxax);
% Efield = (Emax/2)*(cos(cosax)+1.01);
% Efield = [(zeros(1,npts-1-length(cosax))), Efield];
Efield = exp(100*vxax);
Efield = Efield./max(Efield);
Efield = Emax*Efield;
Efield = Efield.^2;
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



