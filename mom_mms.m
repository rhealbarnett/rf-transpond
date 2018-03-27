%----------------------------------------%
% Method of manufactured solutions       %
% mom eqn                                %
% dvx/dt + dvx/dx = 0                    %
% rlbarnett c3149416 230318              %
%----------------------------------------%

%------
% temporal domain %
%------
tmin = 0;

% max number of iterations -- keep constant
nmax = 2000;

% temporal step
dt = 1.0/nmax;

% time axis
tax = linspace(tmin,nmax*dt,nmax);

%------
% spatial domain %
%------
xmin = 0;
xmax = 2.0*pi;

% spatial step and axis
% dx = logspace(-2,0,10);
dx = 1.;

% coefficients for source term s = a + sin(x + ct)
a = 0.0;
c = -1.0;

%%
%-------------------------------------------------------------------------%
% Manufactured solution U(t,x) = sin(x-t)                                 %
% dU/dt = -cos(x-t); d/dx(U^2/2) = -(1/2)sin(2(t-x))                      %
%-------------------------------------------------------------------------%

% diffusion term multiplier
nu = 0.05;

% number of dx iterations (ie times to perform dx=dx/2
iter = 8;
l_sec = zeros(1,iter);
l_inf = zeros(1,iter);

ax = 0;

for kk=1:iter
    
    npts = round(1 + (xmax - xmin)/dx);
    xax = linspace(xmin,xmax,npts);
    
    % initial conditions
    source = zeros(1,npts);
    vx_new = a + sin(xax);
    
    for ii=1:21
        
        % exact solution
        ex_sol = a + sin(xax + c*dt*ii);
        
        vx = vx_new;

        for jj=2:npts-1

            source(1,jj) = c*cos(c*ii*dt + xax(jj)) + (a + sin(xax(jj) + ii*dt*c))*cos(xax(jj) + c*ii*dt) + nu*sin(xax(jj) + c*ii*dt); 

            vx_new(1,jj) = vx(1,jj) - (dt/(2.0*dx))*vx(1,jj)*(vx(1,jj+1) - vx(1,jj-1)) +...
                (dt/(dx^2))*nu*(vx(1,jj+1) - 2.0*vx(1,jj) + vx(1,jj-1)) + dt*source(1,jj);
        end

        % boundary conditions
        source(1,end) = c*cos(xmax + c*dt*ii) + (a + sin(xmax + ii*dt*c))*cos(xmax + c*ii*dt) + nu*sin(xmax + c*ii*dt);
%         vx(1,end) = vx(1,end) - (dt/(2.0*dx))*vx(1,end)*(vx(1,2) - vx(1,end-1)) +...
%                 (dt/(dx^2))*nu*(vx(1,2) - 2.0*vx(1,end) + vx(1,end-1)) + dt*source(1,end);
%         vx(1,1) = vx_new(1,end);
        
        vx(1,end) = a + sin(xmax + c*dt*ii) + source(1,end);
        vx_new(1,1) = vx_new(1,end);

        l_inf(1,kk) = max(abs(ex_sol - vx_new));
        l_sec(1,kk) = rms(ex_sol - vx_new);
       

    end
    
    ax(1,kk) = dx;
    dx = dx/2.0;
    
end


%%
%------ 
% plot solution %
%------

figure(1)
plot(xax, vx, '*k')
hold on
plot(xax, ex_sol)
xlabel('Location')
ylabel('Amplitude')
legend('vx', 'exact solution')
hold off

figure(2)
loglog(ax,l_sec, '-*r') 
hold on
loglog(ax,l_inf, '-ob')
loglog(ax, ax.^2, '--k')
% xlim([1e-2 1e0])
xlabel('dx')
ylabel('Error')
legend('L$_2$', 'L$_{\inf}$', 'dx$^2$','Location','northwest')
hold off
