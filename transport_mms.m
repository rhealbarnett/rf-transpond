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
% npts = 4096;
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
    s = linspace(smin,smax,npts-1); 

    % calculate the x values from xi
    x = xmax*(s.^(1/ro));
%     
%     figure(1);
%     plot(s,x,'LineWidth',2);
%     hold on
%     
%     for ii=1:round(npts/50):npts
%         
%         figure(1)
%         hold on
%         plot([s(ii),s(ii)],[0.0, x(ii)]);
%         plot([0.0,s(ii)],[x(ii), x(ii)]);
%     
%     end
%     
%     xlabel('\xi','Fontsize',16) 
%     ylabel('x','Fontsize',16)
%     set(gca,'Box','on'); 
%     set(gca,'Fontsize',16, 'LineWidth',2)
%     plot([smin,smax],[0,xmax],'r')
%     % text(0.05,L-0.1,['A=',num2str(A),...
%     % ' x=',num2str(xc)],'Fontsize',18)
% 
%     hold off;

    vxax = x;
    vdx = (vxax(2:end) - vxax(1:end-1));
    vdx = fliplr(vdx);
    vxax(1) = xmin;
    for ii=2:(npts/2)-1
        vxax(1,ii) = vxax(1,ii-1) + vdx(1,ii-1);
    end

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

% decay_const = 300.0;

mms_mult = 100000.0;
u0 = 1.0;
% epsilon = 1.0;
epsilon = 0.001;
om = 1.0e6;

%-- initial density profile
n_new = u0*(sin(mms_mult*nxax.^2) + epsilon);
n_avg = interp1(nxax,n_new,vxax);
n_init = n_new;
n_source = zeros(1,npts);

%-- initial velocity
% vx_new = u0*(sin(mms_mult*vxax.^2) + epsilon);
vx_new = exp(-mms_mult*vxax)*(sin(0) + epsilon);
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
% dt = 5.9037e-10;
% dt = 1.4768e-10;
tmin = 0;
tmax = 1.0e-7;
% tmax = 20*dt;
nmax = round(tmax/dt);
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
