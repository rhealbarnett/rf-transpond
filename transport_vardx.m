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
m = 4.0*amu;
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
npts = 2048;
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

    npts = length(vxax) + 1;
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
% n_new = (Nmin*(fliplr(nxax)/max(nxax)) + Nmin);
n_new = Nmax*(1.0 - (5.0e-19)*exp((nxax(end/2 + 1:end))*10));
n_new = [fliplr(n_new), n_new];
n_avg = interp1(nxax,n_new,vxax);
n_init = n_new;

%-- initial velocity
vx_mult = log(cs)/(xmax - xmin);
vx_const = -exp(vx_mult*xmin);
% vx_new = cs*(1.0 + (1.0e-40)*exp((vxax(end/2 + 1:end))*10));
% vx_new = [fliplr(vx_new),vx_new];
% vx_new = u0*(sin(mms_mult*vxax.^2) + epsilon);
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
% 
% rate coefficient (constant)
rate_coeff = (1.0e-14);
% approx size of non-zero portion of neutral profile (1/4 domain)
decay_loc = xmax - 0.0001*xmax;
% decay_loc = xmax - 0.2*xmax;
a = find(nxax >= decay_loc);
decay_index = npts - a(1);

% decay_index = round((npts)/4);
% calculate shape of neutral profile
cosax = linspace(pi,2*pi,decay_index);
% max neutral value (at wall)
neut_max = (1.0e18);

% initialise and fill neutral density array
n_neut = zeros(1,npts);
n_neut(end-decay_index+1:end) = neut_max*((cos(cosax) + 1)/2);
% n_neut(1:end-decay_index+1) = n_neut(end-decay_index+2);
n_neut(end/2 + 1:end-decay_index+1) = neut_max*((cos(pi) + 1)/2);
n_neut(1,1:end/2) = fliplr(n_neut(end/2 + 1:end));
n_neut = interp1(linspace(xmin-0.5*dx,xmax+0.5*dx,npts),n_neut,nxax,'linear');

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
lflux = n_avg(1)*vx_new(1);
% calculate the constant multiplier to match density out = in
ns_mult = (rflux-lflux)/source_int;
% multiply n0(x)n(x,t) by the constant calculated in previous step
n_source = (n_source*ns_mult)*1.0e-2;
nv_source = source_avg*ns_mult;
n_source(1,1) = 0.0; n_source(1,end) = 0.0;

equib = load('../../lapd_numdata/matlab/equibhe_8m.mat');
equib_n_source = equib.n_source;
equib_nxax = equib.nxax;
n_source = interp1(equib_nxax,equib_n_source,nxax,'linear');

% n_source = equib.n_source;

fprintf('Initial outward flux at RH boundary %d\n',rflux)
fprintf('Initial density source integral %d\n',source_int)
fprintf('Density source integral after normalisation %e\n',trapz(nxax,n_source))
fprintf('Initial total number of particles %e\n',trapz(nxax,n_new))
fprintf(['Difference between outward flux at boundaries and integral ' ...
    ' of density source %e\n'], (rflux-lflux) - trapz(vxax,nv_source))    


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
tmax = 1.0e-2;
cfl_fact = 0.99;

% if ((cfl_fact*(min(ndx)^2)/(2.0*nu))<(cfl_fact*min(ndx)/max(abs(vx_new))))
%     dt = cfl_fact*(min(ndx)^2)/(2.0*nu);
% elseif (cfl_fact*min(ndx)^2/(2.0*nu))>(cfl_fact*min(ndx)/max(abs(vx_new)))
%     dt = cfl_fact*min(ndx)/max(abs(vx_new));
% else
%     dt = cfl_fact*min(ndx)/cs;
% end

dt = cfl_fact*min(ndx)/max(abs(vx_new));
dt = 6.0*dt;
% dt = (2.0*pi/om)*0.01;
nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);
% mult = 1.0./dx;

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


