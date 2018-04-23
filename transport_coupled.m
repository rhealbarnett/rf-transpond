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
cs = sqrt((Te + Ti)*e/m);
nu = 1.0;

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 10/cs;
npts = 64;
dx = (xmax - xmin)/(npts - 1);
xax = linspace(xmin,xmax,npts);

%------
% temporal domain %
%------
tmin = 0;

% calculate dt based on vacuum speed of light for now
dt = 8.0e-10;
nmax = 1.0e7;

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
%-------------------------------------------------------------------------%

n_new = normpdf(xax,xmax/2,xmax/30);
n_new = n_new/max(n_new);
n_new = 1.0e3*n_new;
n_new = n_new + 1;
n_new = n_new*1.0e16;
dnx = gradient(n_new,xax);

vx_new = zeros(1,npts);
nA = zeros(npts,npts);
vxA = zeros(npts,npts);
vx_source = zeros(npts,1);

nA(1,1) = 1.0;
nA(end,end) = 1.0;
vxA(1,1) = 1.0;
vxA(end,end) = 1.0;

%%

figure(1)
semilogy(xax,n_new,'DisplayName','time = 0s')
hold on

%-------------------------------------------------------------------------%
% Solve                                                                   %
%-------------------------------------------------------------------------%

for ii=1:nmax

%     if (0.99*(dx^2)/(2.0*nu))<(0.99*dx/max(vx_new))
%         dt = 0.99*(dx^2)/(2.0*nu);
%         continue
%     elseif (0.99*(dx^2)/(2.0*nu))>(0.99*dx/max(vx_new))
%         dt = 0.99*dx/max(vx_new);
%         continue
%     end
    
    vx = vx_new;
    n = n_new;

    for jj=2:npts-1
        
        nA(jj,jj) = 1.0 + (dt/(2.0*dx))*(vx(1,jj-1) - vx(1,jj+1));
        nA(jj,jj-1) = (dt/(2.0*dx))*vx(1,jj);
        nA(jj,jj+1) = -(dt/(2.0*dx))*vx(1,jj);
        
        vxA(jj,jj) = (2.0*nu*dt)/(dx^2) - 1;
        vxA(jj,jj-1) = -(vx(1,jj-1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        vxA(jj,jj+1) = (vx(1,jj+1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        
        vx_source(jj,1) = -((Te + Ti)*e)/(m*n(1,jj)*2.0*dx)*(n(1,jj+1) - n(1,jj-1)); 

    end
    
    nA = sparse(nA);
    vxA = sparse(vxA);
    
    n_new = nA*n';
    vx_new = vx_source - vxA*vx';
    
    n_new = n_new';
    vx_new = vx_new';
    
    vx_new(1,1) = 0;
    vx_new(1,npts) = 0;
    n_new(1,1) = n_new(1,2);
    n_new(1,end) = n_new(1,end-1);
    
    nan_check = isnan(vx_new);
    
    if sum(nan_check) ~= 0
        fprintf('unstable, ii=%d\n',ii)
        break
    elseif dt*max(vx_new)/dx >= 1.0 || dt*2*nu/dx^2 >= 1.0
        fprintf('CFL condition violated')
        break
    else
        continue
    end 
    
    
%     figure(1)
%     plot(xax,n_new,'DisplayName',['time = ' num2str(ii*dt) ' s'])
%     hold on
    
end

%%

semilogy(xax,n_new,'--r','DisplayName',['time = ' num2str(nmax*dt) ' s'])
legend('show')
hold off

figure(2)
semilogy(xax,n_new)
xlabel('Position ($m$)','Fontsize',16)
ylabel('log$_{10}|$Density$|$','Fontsize',16)

figure(3)
plot(xax,vx_new)
xlabel('Position ($m$)','Fontsize',16)
ylabel('Velocity ($ms^{-1}$)','Fontsize',16)




