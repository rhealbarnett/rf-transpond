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
% dx = 0.0013378;
% xax = xmin:dx:xmax;
% npts = length(xax);

%--
% temporal domain
tmin = 0;
dt = 0.99*dx/c0;
tmax = 1.0e7*dt;
nmax = tmax/dt;
tax = linspace(tmin,tmax,nmax);
tol = 1.0e-5;
% dt = 0.99*dx;
% nmax = 756;

%--
% constants
e = 1.6022e-19;
m = 1.67e-27;

%--
% parameters
Te = 5.0;
Ti = 5.0;
T = Te + Ti;
cs = sqrt((Te + Ti)*e/m);

%%

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
n(1,:) = 1.0e19;
% n(1) = 1.0e19;

for ii=1:nmax
    n_old = n;
    for jj=2:npts-1
        n(1,jj) = n(1,jj) - vx(1,jj)*(dt/dx)*(n(1,jj) - n(1,jj-1)) - n(1,jj)*(dt)*dvx(1,jj);%/dx)*(vx(1,jj) - vx(1,jj-1));
    end
    n(1,1) = 1.0e19;
    n(1,end) = n(1,end-1);
%     if abs(rms(n_old)-rms(n))<tol
%         break
%         ii
%     else
%         continue
%     end
end


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
Nmax = 19;
Nmin = 16;
m = (Nmax - Nmin) ./ (xmax - xmin);
n = 10.^(-m*xax + Nmax);
% n = 1.0e19*ones(1,npts);
n = 0.5e19*(xax)+1.05e19;
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

vx = (cs/2)*ones(1,npts);
% vx = zeros(1,npts);
% vx(1,1) = 0;
% vx(1,end) = cs;
% vx = sin(xax);
% a = find(xax<=0.5);
% vx(a(end):end) = 1.0;
vx_new = zeros(1,npts);

figure(5)
plot(xax,vx,'--k')
hold on

for ii=1:nmax
    for jj=2:npts-1
        vx(1,jj) = (1./2.)*(vx(1,jj-1) + vx(1,jj+1)) - (dt/(2.0*dx))*((vx(1,jj+1)^2)/2. - (vx(1,jj-1)^2)/2.) - ((T*e)/(n(1,jj)*m))*(dt/dx)*...
            (n(1,jj) - n(1,jj-1));
    end
    vx(1,1) = 0.0;
    vx(1,end) = cs;
%     vx_new(1,1) = vx(1,1);
%     vx(1,end) = vx(1,end-1);
%     vx = vx_new;
%     if (ii==100) || (ii==300) || (ii==600) || (ii==756)
%         figure(5)
%         plot(xax,vx_new)
%         hold on
%     if abs(rms(vx_old)-rms(vx))<1.0e-9
%         break
%     elseif isnan(vx(end-1))
%         break
%     else
%         continue
%     end
%     end
end

% figure(5)
% legend('t=0s', 't=0.13s', 't=0.4s', 't=0.8s', 't=1.0s','Location','northwest')
% 
% hold off
    
%--
% % plot solution to compare with comsol solution
figure(4)
plot(xax,vx)
xlabel('Position ($m$)','Fontsize',16)
ylabel('Velocity ($ms^{-1}$)','Fontsize',16)
