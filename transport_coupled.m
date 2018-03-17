%----------------------------------------%
% solve coupled transport equations      %
% cont and mom eqns                      %
% dvz/dt + (1/mn)(Te+Ti)dn/dz = 0        %
% dn/dt + (n)dvz/dz = 0                  %
% d/dt convective derivative             %
% rlbarnett c3149416 150318              %
%----------------------------------------%

c0 = 3.0e8;

%--
% spatial domain
xmin = 0;
xmax = 1;
xax = linspace(xmin,xmax,256);
npts = length(xax);
dx = (xmax - xmin)/(npts-1);

%--
% temporal domain
tmin = 0;
dt = 0.99*dx/c0;
tmax = 1.0e6*dt;
nmax = tmax/dt;
tax = linspace(tmin,tmax,nmax);
tol = 1.0e-5;

%--
% constants
e = 1.6022e-19;
m = 1.67e-27;

%--
% parameters
Te = 10.0e3;
Ti = 5.0e3;
T = Te + Ti;
cs = sqrt((Te + Ti)*e/m);

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
%-------------------------------------------------------------------------%

n = zeros(1,npts);
n(1,:) = 1.0e19;

vx = zeros(1,npts);
vx(1) = 0;
vx(npts) = cs;

%%

%-------------------------------------------------------------------------%
% Solve                                                                   %
%-------------------------------------------------------------------------%

for ii=1:nmax
    n_old = n;
    for jj=2:npts-1
        n(1,jj) = n(1,jj) - vx(1,jj)*(dt/dx)*(n(1,jj) - n(1,jj-1)) - n(1,jj)*(dt/dx)*(vx(1,jj) - vx(1,jj-1));
        vx(1,jj) = (1./2.)*(vx(1,jj-1) + vx(1,jj+1)) - (dt/(2.0*dx))*((vx(1,jj+1)^2)/2. - (vx(1,jj-1)^2)/2.) - ((T*e)/(n(1,jj)*m))*(dt/dx)*...
            (n(1,jj) - n(1,jj-1));
    end
    n(1,npts) = n(1,npts-1);
end

figure(1)
semilogy(xax,n)
xlabel('Position ($m$)','Fontsize',16)
ylabel('log$_{10}|$Density$|$','Fontsize',16)

figure(2)
plot(xax,vx)
xlabel('Position ($m$)','Fontsize',16)
ylabel('Velocity ($ms^{-1}$)','Fontsize',16)




