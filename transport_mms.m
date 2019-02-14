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
variable = 1;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% density solution space will be defined as having npts-2 (xax(2:npts-1))
% npts = 32;
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

    % root order
    ro = 2.0;

    % initialise xi parameter array
    % spacing in the centre currently is 0.5*dx
    smax = 1.0;
    smin = 0.0;
    s = linspace(smin,smax,npts-1); 

    % calculate the x values from xi
    x = xmax*(s.^(1/ro));
%     x = exp(s.^1.5) - smax;
%     x = xmax*x/max(x);

    vxax = x;
    vdx = (vxax(2:end) - vxax(1:end-1));
%     vdx = fliplr(vdx);
%     vxax(1) = xmin;
%     for ii=2:npts-1
%         vxax(1,ii) = vxax(1,ii-1) + vdx(1,ii-1);
%     end

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

% mms_mult = 100000.0;

u0 = 5.0;
n0 = 500.0;

ux = 2.0;
nx = 200.0;

knx = 2000.0;
kux = 2000.0;
lamx = 0.0;

% epsilon = 0.001;
om = 1.0e6;
% om = 2.0e9;
om = 0;

%-- initial density profile
LnBC = 0.0;
RnBC = 0.0;
n_new = n0 + nx*sin(knx*nxax.^2 + 0);
LnBC = n0 + nx*sin(knx*min(nxax)^2 + 0);
RnBC = n0 + nx*sin(knx*max(nxax)^2 + 0);
% n_new = n0 + nx*sin(0)*exp(-lamx*nxax);
% LnBC = n0 + nx*sin(0)*exp(-lamx*min(nxax));
% RnBC = n0 + nx*sin(0)*exp(-lamx*max(nxax));
n_init = n_new;
n_avg = interp1(nxax,n_new,vxax);
% n_source = zeros(1,npts);

%-- initial velocity
vx_new = (u0 + ux*cos(kux*vxax.^2 + 0))*0.01;
LuBC = u0 + ux*cos(kux*min(vxax)^2 + 0);
RuBC = u0 + ux*cos(kux*max(vxax)^2 + 0);
% vx_new = u0 + ux*cos(kux*0)*exp(-lamx*vxax);
% LuBC = u0 + ux*cos(kux*0)*exp(-lamx*min(vxax));
% RuBC = u0 + ux*cos(kux*0)*exp(-lamx*max(vxax));
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

cfl_fact = 0.99;

if ((cfl_fact*(min(ndx)^2)/(2.0*nu))<(cfl_fact*min(ndx)/max(abs(vx_new))))
    dt = cfl_fact*(min(ndx)^2)/(2.0*nu);
elseif (cfl_fact*min(ndx)^2/(2.0*nu))>(cfl_fact*min(ndx)/max(abs(vx_new)))
    dt = cfl_fact*min(ndx)/max(abs(vx_new));
else
    dt = cfl_fact*min(ndx)/cs;
end

dt = 2.0*dt;
dt = 1.4768e-10;
tmin = 0;
dt = 0;
tmax = 1000*dt;
nmax = round(tmax/dt);
nmax = 10000;
tol = 1.0e-12;
tax = linspace(tmin,tmax,nmax);

% for ii=1:nmax
%     vx_new = exp(-decay_const*vxax)*(sin(om*dt*ii) + epsilon);
% %     vx_new = 1.0e4*vx_new/max(vx_new);
%     figure(2)
%     plot(vxax,vx_new)
%     pause(2)
%     hold on
% end

% figure(2)
% hold on
% xlabel('Position (m)','Fontsize',16)
% ylabel('Exact solution','Fontsize',16)
% hold off

Efield = zeros(1,npts-1);
