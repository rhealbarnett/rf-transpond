%-------------------------------------------------------%
% MMS file for transport equations                      %
% rlbarnett c3149416 220119                             %
%-------------------------------------------------------%

%------
% constants %
%------
const = constants();
m = const.mp;
e = const.e;

%------
% parameters %
%------
Te = 5.0;
Ti = 10.0;
cs = sqrt((Te + Ti)*e/m);
nu = 1.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 0.1;

%------
% turn variable grid on (1) or off (0)
%------
variable = 0;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% density solution space will be defined as having npts-2 (xax(2:npts-1))
% npts = 4096;
npts = 4096;
dx = (xmax - xmin)/(npts - 1);
nxax = linspace(xmin-0.5*dx,xmax+0.5*dx,npts);
vxax = linspace(xmin,xmax,npts-1);
ndx = dx*ones(1,npts-1);
vdx = dx*ones(1,npts-2);


%%
%----------------------------------------------------------------------%
% non uniform grid calculation
% STILL TESTING
%----------------------------------------------------------------------%

if variable
    
    clear vxax nxax

    % root order
    ro = 2.0;

    % initialise xi parameter array
    % spacing in the centre currently is 0.5*dx
    smax = 1.0;
    smin = 0.0;
%     s = smin:20.0*dx:smax;
    s = linspace(smin,smax,(npts/2)-1); 

    % calculate the x values from xi
    x = xmax*(s.^(1/ro));

    vxax = x;
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

mms_mult = 100.00;
u0 = 1.0;
epsilon = 0.001;
om = 1.0e6;

%-- initial density profile
n_new = u0*(sin(mms_mult*nxax.^2) + epsilon);
n_avg = interp1(nxax,n_new,vxax);
n_init = n_new;

%-- initial velocity
vx_new = u0*(sin(mms_mult*vxax.^2) + epsilon);
vx_init = vx_new;


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

tmin = 0;
tmax = 1.0e-5;
cfl_fact = 0.99;

if ((cfl_fact*(min(ndx)^2)/(2.0*nu))<(cfl_fact*min(ndx)/max(abs(vx_new))))
    dt = cfl_fact*(min(ndx)^2)/(2.0*nu);
elseif (cfl_fact*min(ndx)^2/(2.0*nu))>(cfl_fact*min(ndx)/max(abs(vx_new)))
    dt = cfl_fact*min(ndx)/max(abs(vx_new));
else
    dt = cfl_fact*min(ndx)/cs;
end

dt = 2.0*dt;
nmax = round(tmax/dt);
tax = linspace(tmin,tmax,nmax);
