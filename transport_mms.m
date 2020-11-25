%-------------------------------------------------------%
% MMS file for transport equations                      %
% rlbarnett c3149416 220119                             %
%-------------------------------------------------------%

%--
% Constants 
const = constants();
m = const.mp;
e = const.e;

%--
% Parameters 
Te = 5.0;
Ti = 10.0;
cs = sqrt((Te + Ti)*e/m);
nu = 0.7;

%--
% Spatial domain
xmin = -0.3;
xmax = 0.4;

%--
% Turn variable grid on (1) or off (0)
variable = 1;

% Define uniform grids. For a staggered grid, n has npts and v has npts-1.
% The number of points (npts) should be set in the run_mms.m function.
dx = (xmax - xmin)/(npts);
nxax = linspace(xmin-0.5*dx,xmax+0.5*dx,npts);
vxax = linspace(xmin,xmax,npts-1);
ndx = dx*ones(1,npts-1);
vdx = dx*ones(1,npts-2);


%%
%----------------------------------------------------------------------%
% Non-uniform grid calculation
%----------------------------------------------------------------------%

if variable
    
    clear vxax nxax
    
    %--
    % Location of sign change
    xc = (xmax - xmin)/2.0;
    
    %--
    % 'Strength' of grid refinement.
    % Sign also indicates whether refinement is in the centre or at the
    % boundaries.
    A = -1.0;

    %--
    % Set up the unit spaced parameter, xi, that the grid is a function of.
    smax = 1.0;
    smin = 0.0;
    s = linspace(smin,smax,npts-1); 

    %--
    % Calculate the x values from xi.
    x = xc*(1.0 - tanh(A*(1.0 - 2.0*s))./tanh(A));

    %--
    % Set vxax based on the variable x grid, shifted by xmin.
    vxax = x - abs(xmin);
    vdx = (vxax(2:end) - vxax(1:end-1));
    
    %--
    % Reinitialise nxax and calculate based on vxax. 
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

u0 = 0.01;
n0 = 0.01;

ux = 1.0;
nx = 1.0;

knx = 1.0;
kux = 5.0;

%-- initial density profile
if SS
    n_new = 0.01*(n0 + nx*sin(knx*nxax.^2 + 0));
    vx_new = 0.01*(u0 + ux*cos(kux*vxax.^2 + 0));
    
    om = 0.0;
elseif TD
    n_new = (n0 + nx*sin(knx*nxax.^2 + 0));
    vx_new = (u0 + ux*cos(kux*vxax.^2 + 0));
    
    om = 1.0e5;
end

LnBC = n0 + nx*sin(knx*min(nxax)^2 + 0);
RnBC = n0 + nx*sin(knx*max(nxax)^2 + 0);
n_init = n_new;
n_avg = interp1(nxax,n_new,vxax);
ex_soln = n0 + nx*sin((nxax).^2 + 0);


%-- initial velocity
LuBC = u0 + ux*cos(kux*min(vxax)^2 + 0);
RuBC = u0 + ux*cos(kux*max(vxax)^2 + 0);
vx_init = vx_new;
ex_solu = u0 + ux*cos(kux*vxax.^2 + 0);

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

rf_ex = zeros(1,npts-1);
rf_ey = zeros(1,npts-1);
rf_ez = zeros(1,npts-1);
