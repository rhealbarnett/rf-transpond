%--------------------------------------------------------------------------------------------------------------%
% Input file for transport 1.5D                                                                            %
% nabla = (kapx, kapy, dz)
% rlbarnett c3149416 190723     
%--------------------------------------------------------------------------------------------------------------%
%%

%-------
% call file containing physical constants. 
%-------
const = constants();
mp = const.mp;
e = const.e;
mu0 = const.mu0;

%-------
% define parameters. 
%-------
Te = 5.0;
Ti = 1.0;
cs = sqrt((Te + Ti)*e/mp);
B0 = 1.0;

eta_para = 1.0;

%-- 
% D_perp will be used to calculate eta_perp once the density (n) 
% has been defined. 
D_perp = 1.0;

%-- 
% Setting the perp derivatives as zero to begin, for comparison with the 
% parallel only case. 
kap_x = 0.1;
kap_y = 0.1;

%%
%-------
% Initialise the grid parameters.
%-------

npts = 512;
xmin = -4.0;
xmax = 4.0;

variable = 1;

%--------
% Calculate variable grid OR keep uniform grid.
%--------

if variable
    %-- 
    % Find location of sign change.
    xc = (xmax - xmin)/2.0;

    %--
    % Strength of grid refinement. Sign indicates whether refinement is in 
    % the centre (positive) or at the boundaries (negative).
    A = -5.0;

    %--
    % set up the unit spaced parameter, s, that the grid is a function of.
    smax = 1.0;
    smin = 0.0;
    s = linspace(smin,smax,npts-1); 

    %--
    % Calculate the grid values from s. 
    x = xc*(1.0 - tanh(A*(1.0 - 2.0*s))./tanh(A));

    %--
    % Calculate the velocity grid & grid spacing from the calculated grid above. 
    vxax = x - xmax;
    vdx = (vxax(2:end) - vxax(1:end-1));

    %--
    % Initialise density grid.
    nxax = zeros(1,npts);

    %--
    % Fill first (last) index with value calculated from first (last) velocity 
    % grid value. 
    nxax(1) = vxax(1) - 0.5*vdx(1);
    nxax(end) = vxax(end) + 0.5*vdx(end);

    %--
    % Calculate the density grid as halfway between velocity grid points. Also
    % calculate density grid spacing. 
    nxax(2:end-1) = (vxax(2:end) + vxax(1:end-1))/2.0;
    ndx = nxax(2:end) - nxax(1:end-1);
    
elseif ~variable
    
    %-- 
    % Calculate grid spacing.
    dx = (xmax - xmin)/(npts - 1);
    
    %-- 
    % Define density and velocity grids.
    nxax = linspace(xmin-0.5*dx,xmax+0.5*dx,npts);
    vxax = linspace(xmin,xmax,npts-1);
    
    %-- 
    % Calculate uniform spaced grid. Needs to have this form to maintain
    % consistency with the transport solve itself.
    ndx = dx*ones(1,npts-1);
    vdx = dx*ones(1,npts-2);
    
end

%%
%-------
% Calculate perpendicular drift velocity from analytic expression
% (in this case, PF only).
%-------

%--
% Include analytic expression for the electric field. 
mu = 0.0;
sigma = 1.0;
Ex = zeros(1,npts-1);
Ey = (1.0/sqrt(2*pi*sigma^2))*exp(-(vxax - mu).^2/(2.0*sigma^2));
Ey = (Ey ./ (max(Ey)))*300;
Ez = zeros(1,npts-1);
E = [Ex; Ey; Ez];

%--
% Initialise a magnetic field vector. NOTE: this is currently on a uniform
% grid vs a variable grid for the acceleration. Okay for now, as the B
% field is constant. 
B = [zeros(1,npts-1); zeros(1,npts-1); B0*ones(1,npts-1)];

%--
% Calculate v drift from F X B0 
vdrift = cross(E,B);
vdrift_x = vdrift(1,:);
vdrift_y = vdrift(2,:);


%%
%-------
% Initialise density and velocity.
%-------

%-- 
% Max and min density values for the initial profile.
Nmax = 1.0e18;
Nmin = 0.5e18;

%--
% Calculate initial density profile. 
n_new = Nmax*(1.0 - (5.0e-19)*exp((nxax(end/2 + 1:end))*10));
n_new = [fliplr(n_new), n_new];
n_avg = interp1(nxax,n_new,vxax);
n_init = n_new;

%--
% Velocity boundary conditions. Set to +- the sound speed for metal walls. 
LvBC = -cs;
RvBC = cs;

%-- 
% Initial velocity
vx_new = zeros(1,npts-1);
vx_new(1,1) = LvBC;
vx_new(1,end) = RvBC;
vx_init = vx_new;

%%
%-------
% Calculate perpendicular viscosity coefficient using density.
%-------

%--
% Starting with a constant estimated by one density value.
eta_perp = D_perp*mp*n_init;

%%
%-------
% Initialise coefficient matrices.                                         
%-------

%-- Initialise all coefficient matrices as sparse matrices. 
nA = sparse(npts,npts);
nI = sparse(eye(npts,npts));
vx_pos = sparse(npts-1,npts-1);
vx_neg = sparse(npts-1,npts-1);
vx_diff = sparse(npts-1,npts-1);
vx_I = sparse(eye(npts-1,npts-1));

%%
%-------
% Calculate time step                                                                  
%-------

%-- 
% Explicit convective dt based on CFL conditions, check during loop if violated
tmin = 0.0;
tmax = 1.0e-7;

%--
% CFL condition multiplier. Keep close to unity, or dt will be very small.  
cfl_fact = 0.99;

%--
% Calculate using the sound speed as that should be the maximum velocity. 
% (Is that true? Or do we need to consider the Alfven velocity?).
dt = cfl_fact*min(ndx)/cs;

%--
% Number of iterations (nmax) and time axis. 
nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);

%%
%-------
% Calculate density source term.
%-------

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

































