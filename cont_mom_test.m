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
xmax = 2.0;
dx = 1.125e-4;
% npts = 128;
npts = round((xmax - xmin)/(dx) + 1);
% dx = (xmax - xmin)/(npts-1);
xax = linspace(xmin,xmax,npts);

%--
% constants
e = 1.6022e-19;
m = 1.67e-27;

%--
% parameters
Te = 10.0;
Ti = 5.0;
T = Te + Ti;
% cs = sqrt((Te + Ti)*e/m);
cs = 1.0e4;
nu = 1.0;

%--
% temporal domain
tmin = 0;
% dt = 0.99*dx/c0;
% dt = 8.0e-6;
dt = 0.99*(dx^2)/(2.0*nu);
tmax = 1.0e5*dt;
nmax = tmax/dt;
tol = 1.0e-14;

%%

% %-------------------------------------------------------------------------%
% % Solve continuity equation with analytic velocity expression             %
% %-------------------------------------------------------------------------%
% 
% %--
% % velocity
% vx = cs*xax;
% dvx = gradient(vx,xax);
% 
% %--
% % plot velocity to compare against comsol output
% figure(1)
% plot(xax,vx)
% hold on
% plot(xax,dvx)
% legend('vx', 'dvx')
% xlabel('Position ($m$)')
% hold off
% 
% %--
% % up wind differencing scheme
% 
% n_new = 1.0e19*ones(1,npts);
% 
% %--
% % initialise matrix
% coeff_mat = zeros(npts,npts);
% coeff_mat(1,1) = 1.0;
% coeff_mat(end,end) = 1.0;
% 
% for ii=1:nmax
%     
%     n = n_new;
%     
%     for jj=2:npts-1
%         
%         coeff_mat(jj,jj) = (dt/(2.0*dx))*(vx(1,jj-1) - vx(1,jj+1));
%         coeff_mat(jj,jj-1) = (1.0/2.0) + (dt/(2.0*dx))*vx(1,jj);
%         coeff_mat(jj,jj+1) = (1.0/2.0) - (dt/(2.0*dx))*vx(1,jj);
%         
%     end
%     
%     n_new = coeff_mat*n';
%     
%     n_new = n_new';
%     
%     n_new(1,1) = 1.0e19;
%     n_new(1,end) = n_new(1,end-1);
%     
%     if rms(n - n_new)<=tol
%         fprintf('tolerance reached, ii=%d\n',ii)
%         break
%     else
%         continue
%         
%     end
% 
% end
% 
% 
% %--
% % plot solution to compare with comsol solution
% figure(2)
% semilogy(xax,n)
% xlabel('Position ($m$)','Fontsize',16)
% ylabel('log$_{10}|$Density$|$','Fontsize',16)

%%
%-------------------------------------------------------------------------%
% Solve momentum equation with analytic density expression                %
%-------------------------------------------------------------------------%

clear vx dvx n

%--
% density
% Nmax = 1.0e19;
% n = nmax*ones(1,npts);
n = normpdf(xax,1,0.05);
n = n/max(n);
n = 1.0e3*n;
n = n + 1;
n = n*1.0e16;
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
% initialise velocity coefficient matrix and rhs
vxA = zeros(npts,npts);
vx_source = zeros(npts,1);

%--
% initial value for velocity
vx_new = (cs/4.0)*xax + cs/2.0;
% vx_new = cs*xax;

for ii=1:nmax
    
    vx = vx_new;

    for jj=2:npts-1
        
        vxA(jj,jj) = (2.0*nu*dt)/(dx^2) - 1.0;
        vxA(jj,jj-1) = -(vx(1,jj-1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        vxA(jj,jj+1) = (vx(1,jj+1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        
        vx_source(jj,1) = 0.0;%-((Te + Ti)*e)/(m*n(1,jj)*2.0*dx)*(n(1,jj+1) - n(1,jj-1)); 

    end
    
    vxA = sparse(vxA);
    
    vx_new = vx_source - vxA*vx';
  
    vx_new(1,1) = cs/2;
    vx_new(end,1) = cs;
   
    vx_new = vx_new';
    
    nan_check = isnan(vx_new);
    
    if sum(nan_check) ~= 0
        fprintf('unstable, ii=%d\n',ii)
        break
    else
        continue
    end
    
    if (rms(vx - vx_new) <= tol)
        fprintf('tolerance reached, ii=%d\n',ii)
        break
    else
        continue
        
    end
end

%%
%--

figure(6)
plot(xax,vx_source) 
hold on
plot(xax,vx_new)
legend('source','vx')
hold off
