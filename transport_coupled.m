%----------------------------------------%
% solve coupled transport equations      %
% cont and mom eqns                      %
% dvz/dt + (1/mn)(Te+Ti)dn/dz = 0        %
% dn/dt + (n)dvz/dz = 0                  %
% d/dt convective derivative             %
% rlbarnett c3149416 150318              %
%----------------------------------------%

%------
% spatial domain %
%------
xmin = 0.0;
xmax = 1.0;
npts = 128;
dx = (xmax - xmin)/(npts - 1);
xax = linspace(xmin,xmax,npts);

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
% temporal domain %
%------
tmin = 0;

% calculate dt based on vacuum speed of light for now
dt = 0.99*dx/c0;
nmax = 1.0e5;
tol = 1.0e-14;

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
%-------------------------------------------------------------------------%

n_new = 1.0e19*ones(1,npts);

vx_new = cs*xax;
nA = zeros(npts,npts);
vxA = zeros(npts,npts);
vx_source = zeros(npts,1);

nA(1,1) = 1.0;
nA(end,end) = 1.0;
vxA(1,1) = 1.0;
vxA(end,end) = 1.0;

%%

%-------------------------------------------------------------------------%
% Solve                                                                   %
%-------------------------------------------------------------------------%

for ii=1:nmax
    
    vx = vx_new;
    n = n_new;

    for jj=2:npts-1
        
        nA(jj,jj) = 1.0 + (dt/(2.0*dx))*(vx(1,jj-1) - vx(1,jj+1));
        nA(jj,jj-1) = (dt/(2.0*dx))*vx(1,jj);
        nA(jj,jj+1) = -(dt/(2.0*dx))*vx(1,jj);
        
        vxA(jj,jj) = (2.0*nu*dt)/(dx^2) - 1;
        vxA(jj,jj-1) = -(vx(1,jj-1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        vxA(jj,jj+1) = (vx(1,jj+1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        
        vx_source(jj,1) = 0;%-((Te + Ti)*e)/(m*n(1,jj)*2.0*dx)*(n(1,jj+1) - n(1,jj-1)); 

    end
    
    nA = sparse(nA);
    vxA = sparse(vxA);
    
    n_new = nA*n';
    vx_new = vx_source - vxA*vx';
    
    n_new = n_new';
    vx_new = vx_new';
    
    vx_new(1,1) = 0;
    vx_new(1,npts) = cs;
    n_new(1,1) = 1.0e19;
    n_new(1,end) = n_new(1,end-1);
    
    if (rms(n - n_new)<=tol) || (rms(vx - vx_new)<=tol)
        fprintf('tolerance reached, ii=%d\n',ii)
        break
    else
        continue
        
    end
end

figure(1)
semilogy(xax,n)
xlabel('Position ($m$)','Fontsize',16)
ylabel('log$_{10}|$Density$|$','Fontsize',16)

figure(2)
plot(xax,vx)
xlabel('Position ($m$)','Fontsize',16)
ylabel('Velocity ($ms^{-1}$)','Fontsize',16)




