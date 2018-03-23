%----------------------------------------%
% test comsol solution 'decoupled'       %
% cont and mom eqns                      %
% dvz/dt + (1/mn)(Te+Ti)dn/dz = 0        %
% dn/dt + (n)dvz/dz = 0                  %
% rlbarnett c3149416 070318              %
%----------------------------------------%

c0 = 3.0e8;

%--
% temporal domain
tmin = 0;
% dt = 0.99*(dx);
% dt = logspace(-5,-11,7);
dt = 1.5e-13;
tmax = 1.0e4*dt;
% tmax = 2.5e-4;
nmax = tmax/dt;
% tax = linspace(tmin,tmax,nmax);
tol = 1.0e-5;
% dt = 0.99*dx;
% nmax = 756;

%--
% spatial domain
xmin = 0;
xmax = 1;
dx = logspace(-1,-5,5);
% npts = (xmax-xmin)./dx;
% xax = linspace(xmin,xmax,npts);
% npts = length(xax);
% dx = (xmax - xmin)/(npts-1);
% dx = 0.0013378;
% xax = xmin:dx:xmax;
% npts = length(xax);

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
% cs = 10.0;



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

%%
%-------------------------------------------------------------------------%
% Solve momentum equation with analytic density expression                %
%-------------------------------------------------------------------------%

clear vx dvx n
% 
% %--
% % density
% Nmax = 19;
% Nmin = 16;
% m = (Nmax - Nmin) ./ (xmax - xmin);
% n = 10.^(-m*xax + Nmax);
% n = 1.0e19*ones(1,npts);
% % n = 0.5e19*(xax)+1.05e19;
% dnx = gradient(n,xax);
% 
% %--
% % plot density to compare against comsol output
% figure(3)
% plot(xax,n)
% hold on
% plot(xax,dnx)
% legend('n', 'dnx')
% xlabel('Position ($m$)')
% hold off

%--
% up wind differencing scheme
% CONSERVATIVE FORM ie (all partials) dv/dt + d(v^2/2)/dx = ...

% vx = (cs)*ones(1,npts);
% vx = zeros(1,npts);
% vx(1,1) = 0;
% vx(1,end) = cs;
% vx = sin(xax);
% % a = find(xax<=0.5);
% % vx(a(end):end) = 1.0;
% vx_new = zeros(1,npts);
% source = zeros(1,npts);
nu = 1.0e-5;

% figure(4)
% plot(xax,vx,'--k')
% hold on


for kk=length(dx)
    kk
    for ii=1:nmax
        npts = round((xmax-xmin)/dx(kk));
        xax = linspace(xmin,xmax,npts);
        vx = sin(xax);
        source = zeros(1,npts);
        l_inf = zeros(1,npts);
        l_sec = zeros(1,npts);
%         tax = linspace(tmin,tmax,round(nmax(kk)));
        for jj=2:npts-2
            vx(1,jj) = (1./2.)*(vx(1,jj-1) + vx(1,jj+1)) - vx(1,jj)*(dt/(2.0*dx(kk)))*((vx(1,jj+1)) - (vx(1,jj-1))) + ...
                (dt*nu/(dx(kk)^2))*(vx(1,jj) - 2.0*vx(1,jj+1) + vx(1,jj+2));
%             source(1,jj) = - cos(ii*dt(kk) - xax(jj)) - (1.0/2.0)*sin(2.0*(ii*dt(kk) - xax(jj))); 
            source(1,jj) = - cos(ii*dt - jj*dx(kk)) + sin(ii*dt - jj*dx(kk))*cos(jj*dx(kk) - ii*dt) + nu*sin(jj*dx(kk) - ii*dt);
%             source(1,jj) = - exp(xax(jj) - ii*dt(kk)) + exp(2.0*(xax(jj) - ii*dt(kk)));
        end
        vx(1,1) = sin(-ii*dt);
        vx(1,end) = sin(npts*dx(kk) - ii*dt);
        source(1,1) = sin(-ii*dt);
        source(1,end) = sin(npts*dx(kk) - ii*dt);
        l_inf(1,ii) = abs(max(vx - source));
        l_sec(1,ii) = rms(vx - source);
        
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
end
% figure(5)
% legend('t=0s', 't=0.13s', 't=0.4s', 't=0.8s', 't=1.0s','Location','northwest')
% 
% hold off
    
%--
% % plot solution to compare with comsol solution
figure(5)
plot(xax,vx)
xlabel('Position ($m$)','Fontsize',16)
ylabel('Velocity ($ms^{-1}$)','Fontsize',16)

figure(6)
plot(xax,source) 
hold on
plot(xax,vx)
legend('source','vx')
