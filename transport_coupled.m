%----------------------------------------%
% solve coupled transport equations      %
% cont and mom eqns                      %
% dvz/dt + (1/mn)(Te+Ti)dn/dz = 0        %
% dn/dt + (n)dvz/dz = 0                  %
% d/dt convective derivative             %
% rlbarnett c3149416 150318              %
%----------------------------------------%

%------
% constants %
%------
e = 1.6022e-19;
c0 = 3.0e8;
m = 1.67e-27;

%------
% parameters %
%------
Te = 10.0;
Ti = 5.0;
T = Te + Ti;
cs = 10*sqrt((Te + Ti)*e/m);
nu = 1.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 10/cs;
npts = 1024;
dx = (xmax - xmin)/(npts - 1);
xax = linspace(xmin,xmax,npts);

%------
% temporal domain %
%------
tmin = 0;

% set dt based on CFL conditions, check during loop if violated
dt = 1.0e-16;
nmax = 1.0e4;
alpha = dt/(dx);

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
%-------------------------------------------------------------------------%

n_new = normpdf(xax(1:end-1),xmax/2,xmax/30);
n_new = n_new/max(n_new);
n_new = 1.0e3*n_new;
n_new = n_new + 1;
n_new = n_new*1.0e16;
dnx = gradient(n_new,xax(1:end-1));

vx_new = zeros(1,npts);
vx_ax = linspace(0,1,npts);
vx_new = -(cs/2)*vx_ax + cs;
% vx_new = cs*ones(1,npts);
% vx_new(1:end/2+1) = -cs;
% vx_new(end/2+1:end) = cs;

nA = zeros(npts-1,npts-1);
vxA = zeros(npts,npts);
vx_source = zeros(npts,1);
phi_plus = zeros(1,npts);
phi_minus = zeros(1,npts);

nA(1,1) = 1.0;
nA(end,end) = 1.0;
vxA(1,1) = 1.0;
vxA(end,end) = 1.0;

%%

figure(1)
plot(xax(1:end-1),n_new,'DisplayName','time = 0s')
ylim([-1.0e19 1.0e19])
% pause(1)
hold on

figure(2)
plot(xax,vx_new,'DisplayName','time = 0s')
hold on

count = 1;

%-------------------------------------------------------------------------%
% Solve                                                                   %
%-------------------------------------------------------------------------%

for ii=1:nmax

%     if (0.99*(dx^2)/(2.0*nu))<(0.99*dx/max(vx_new))
%         dt = 0.99*(dx^2)/(2.0*nu);
%     elseif (0.99*(dx^2)/(2.0*nu))>(0.99*dx/max(vx_new))
%         dt = 0.99*dx/max(vx_new);
%     end
    
    vx = vx_new;
    n = n_new;

    for jj=2:npts-2
        
        if ((vx(1,jj+1)+vx(1,jj))/2)>0
            nA(jj,jj) = 1.0 - alpha*vx(1,jj+1);
            nA(jj,jj-1) = alpha*vx(1,jj);
            nA(npts-1,npts-1) = 1.0 - alpha*vx(1,npts);
            nA(npts-1,npts-2) = alpha*vx(1,npts-1);
        elseif ((vx(1,jj+1)+vx(1,jj))/2)<0
            nA(jj,jj) = 1.0 + alpha*vx(1,jj);
            nA(jj,jj+1) = -alpha*vx(1,jj+1);
            nA(1,1) = 1.0 + alpha*vx(1,1);
            nA(1,2) = -alpha*vx(1,2);
        end        
    end
    
    for jj=2:npts-1
        
        vxA(jj,jj) = (2.0*nu*dt)/(dx^2) - 1;
        vxA(jj,jj-1) = -(vx(1,jj-1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        vxA(jj,jj+1) = (vx(1,jj+1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        
        vx_source(jj,1) = -((Te + Ti)*e)/(m*n(1,jj)*dx)*(n(1,jj) - n(1,jj-1)); 

    end
    
    nA = sparse(nA);
    vxA = sparse(vxA);
    
    n_new = nA*n';
    vx_new = dt*vx_source - vxA*vx';
    
    n_new = n_new';
    vx_new = vx_new';
    
    vx_new(1,1) = cs;
    vx_new(1,end) = cs/2;
    
    nan_check = isnan(vx_new);
    
    if sum(nan_check) ~= 0
        fprintf('unstable, ii=%d\n',ii)
        break
    elseif dt*max(vx_new)/dx >= 1.0 || dt*2*nu/dx^2 >= 1.0
        fprintf('CFL condition violated, ii=%d\n',ii)
        break
    elseif ii==count*1000
        figure(1)
        plot(xax(1:end-1),n_new,'DisplayName',['time = ' num2str(ii*dt) ' s'])
        ylim([-1.0e19 1.0e19]) 
%         legend('show')
        hold on
%         pause(1)
        figure(2)
        plot(xax,vx_new,'DisplayName',['time = ' num2str(ii*dt) ' s'])
        hold on
        count = count + 1;
    end 
    
end

%%

figure(1)
legend('show','Location','northwest')
xlabel('Position (m)')
ylabel('Density m$^{-3}$')
hold off

figure(2)
legend('show','Location','northwest')
xlabel('Position (m)')
ylabel('Velocity ms$^{-1}$')
hold off


figure(3)
plot(xax,vx_new)
xlabel('Position ($m$)','Fontsize',16)
ylabel('Velocity ($ms^{-1}$)','Fontsize',16)




