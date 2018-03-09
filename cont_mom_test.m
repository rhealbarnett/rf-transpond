%----------------------------------------%
% test comsol solution 'decoupled'       %
% cont and mom eqns                      %
% dvz/dt + (1/mn)(Te+Ti)dn/dz = 0        %
% dn/dt + (n)dvz/dz = 0                  %
% rlbarnett c3149416 070318              %
%----------------------------------------%

c0 = 3.0e8;

%--
% spatial domain
xmin = 0;
xmax = 2;
xax = linspace(xmin,xmax,256);
npts = length(xax);
dx = (xmax - xmin)/(npts-1);

%--
% temporal domain
tmin = 0;
dt = 0.99*dx/c0;
tmax = 1.0e7*dt;
nmax = tmax/dt;
tax = linspace(tmin,tmax,nmax);
tol = 1.0e-5;

%--
% constants
e = 1.6022e-19;
m = 1.67e-27;

%--
% parameters
Te = 5.0;
Ti = Te;
T = Te + Ti;
cs = sqrt((Te + Ti)*e/m);

%-------------------------------------------------------------------------%
% Solve continuity equation with analytic velocity expression             %
%-------------------------------------------------------------------------%

%--
% velocity
vx = cs*xax;
dvx = gradient(vx,xax);

%--
% plot velocity to compare against comsol output
figure(1)
plot(xax,vx)
hold on
plot(xax,dvx)
legend('vx', 'dvx')
xlabel('Position ($m$)')
hold off

%--
% up wind differencing scheme

n = zeros(1,npts);
n(1) = 1.0e19;

% for ii=1:nmax
%     n_old = n;
%     for jj=2:npts-1
%         n(1,jj) = n(1,jj) - vx(1,jj)*(dt/dx)*(n(1,jj) - n(1,jj-1)) - n(1,jj)*dt*dvx(1,jj);
%     end
%     n(1,npts) = n(1,npts-1);
%     if abs(rms(n_old)-rms(n))<tol
%         break
%         ii
%     else
%         continue
%     end
% end


%--
% plot solution to compare with comsol solution
figure(2)
semilogy(xax,n)
xlabel('Position ($m$)','Fontsize',16)
ylabel('log$_{10}|$Density$|$','Fontsize',16)

%-------------------------------------------------------------------------%
% Solve momentum equation with analytic density expression                %
%-------------------------------------------------------------------------%
%%
clear vx dvx n

%--
% density
n = 0.5e19*xax+1.0e17;
dnx = gradient(n,xax);

%--
% plot density to compare against comsol output
figure(3)
plot(xax,n)
hold on
plot(xax,dnx)
legend('n', 'dnx')
xlabel('Position ($m$)')
hold off

%--
% up wind differencing scheme
% CONSERVATIVE FORM ie (all partials) dv/dt + d(v^2/2)/dx = ...

vx = zeros(1,npts);
vx(1) = cs;
vx(npts) = 0;

for ii=1:nmax
    vx_old = vx;
    for jj=2:npts-1
        vx(1,jj) = vx(1,jj) - (dt/(2.0*dx))*(vx(1,jj)^2 - vx(1,jj-1)^2) - ((T*e)/(n(1,jj)*m))*dt*dnx(1,jj);
    end
%     if abs(rms(vx_old)-rms(vx))<1.0e-9
%         break
%     elseif isnan(vx(end-1))
%         break
%     else
%         continue
%     end
end