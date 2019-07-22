%--------------------------------------------------------------------------------------------------------------%
% Input file for transport 1.5D                                                                            %
% nabla = (kapx, kapy, dz)
% rlbarnett c3149416 190723     
%--------------------------------------------------------------------------------------------------------------%
%%

%-------
% call file containing physical constants. %
%-------
const = constants();
mp = const.mp;
e = const.e;

%-------
% define parameters. %
%-------
Te = 5.0;
Ti = 1.0;
cs = sqrt((Te + Ti)*e/m);

eta_para = 1.0;

%-- 
% D_perp will be used to calculate eta_perp once the density (n) 
% has been defined. 
D_perp = 1.0;

%-- 
% Setting the perp derivatives as zero to begin, for comparison with the 
% parallel only case. 
kap_x = 0.0;
kap_y = 0.0;

%%
%-------
% Calculate perpendicular drift velocity from analytic expression
% (in this case, PF only).
%-------

%--
% Setting drift velocity as zero to begin, also for comparison with the
% parallel only case. 
vdrift_x = 0.0;
vdrift_y = 0.0;

%%
%-------
% Initialise the grid parameters.
%-------

npts = 512;
xmin = -0.1;
xmax = 0.1;

variable = 0;

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
LuBC = -cs;
RuBC = cs;

%-- 
% Initial velocity
vx_new = zeros(1,npts-1);
vx_new(1,1) = LuBC;
vx_new(1,end) = RuBC;
vx_init = vx_new;



































