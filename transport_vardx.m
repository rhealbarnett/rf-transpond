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
const = constants();
mp = const.mp;
amu = const.amu;
m = mp;
e = const.e;

%------
% parameters %
%------
Te = 5.0;
Ti = 0.5;
cs = sqrt((Te + Ti)*e/m);
nu = 1.0;

%------
% spatial domain %
%------
xmin = -4.0;
xmax = 4.0;

%------
% turn variable grid on (1) or off (0)
%------
variable = 1;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 512;
dx = (xmax - xmin)/(npts - 1);
nxax = linspace(xmin-0.5*dx,xmax+0.5*dx,npts);
vxax = linspace(xmin,xmax,npts-1);
ndx = dx*ones(1,npts-1);
vdx = dx*ones(1,npts-2);


%%
%----------------------------------------------------------------------%
% non uniform grid calculation
%----------------------------------------------------------------------%

if variable
    
    clear vxax nxax

    % location of sign change
    xc = (xmax - xmin)/2.0;
    
    %'strength' of grid refinement.
    % sign also indicates whether refinement is in the centre or at the
    % boundaries
    A = -5.0;

    % set up the unit spaced parameter, xi, that the grid is a function of
    smax = 1.0;
    smin = 0.0;
    s = linspace(smin,smax,npts-1); 

    % calculate the x values from xi
    x = xc*(1.0 - tanh(A*(1.0 - 2.0*s))./tanh(A));

    vxax = x - xmax;
    vdx = (vxax(2:end) - vxax(1:end-1));

%     npts = length(vxax) + 1;
    nxax = zeros(1,npts);

    nxax(1) = vxax(1) - 0.5*vdx(1);
    nxax(end) = vxax(end) + 0.5*vdx(end);
    nxax(2:end-1) = (vxax(2:end) + vxax(1:end-1))/2.0;
    ndx = nxax(2:end) - nxax(1:end-1);
    
end


%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
%-------------------------------------------------------------------------%

% equib = load('equib2.mat');

%-- initial density profile
Nmax = (1.0e18);
Nmin = (0.5e18);
n_new = Nmax*(1.0 - (5.0e-19)*exp((nxax(end/2 + 1:end))*10));
n_new = [fliplr(n_new), n_new];
n_avg = interp1(nxax,n_new,vxax);
n_init = n_new;

%-- initial velocity
LuBC = -cs;
RuBC = cs;
vx_new = zeros(1,npts-1);
vx_new(1,1) = LuBC;
vx_new(1,end) = RuBC;
vx_init = vx_new;


%%
%-------------------------------------------------------------------------%
% CALCULATE DENSITY SOURCE                                                %
%-------------------------------------------------------------------------%

%--
% Rate coefficient (need to double check this value for each ion species
% and temperature... which textbook is this from? Stangeby?).

rate_coeff = (1.0e-14);

%--
% Maximum neutral density value.

neut_max = Nmax;

%--
% Call function to calculate the density source term. 

n_source = density_source(rate_coeff,0.1,nxax,vxax,npts,neut_max,vx_init,n_init); 


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
tmin = 0.0;
tmax = 5.0e-3;
cfl_fact = 0.99;

dt = cfl_fact*min(ndx)/max(abs(vx_new));

nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);


%%

%-----------------------------------------------------------------------------%
% Set up electric field profile & bits for PF                                 %
% Using E||,max = 300Vcm^-1 from J. Myra 2006 Nonlinear interactions paper    %
% exponential decay away from antenna location                                %
%-----------------------------------------------------------------------------%

Emax = 6.0e4;
freq = 80.0e6;
% om = 2.0*pi*freq;
om = 0.0;

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


