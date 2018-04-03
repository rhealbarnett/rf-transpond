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
xmin = -0.1;
xmax = 0.7;

% initial grid size
npts = 11;

% coefficients for source term u0(sin(x^2 + omt) + eps) and viscosity
u0 = 1.0;
eps = 0.001;
nu = 0.7;

%%
%-------------------------------------------------------------------------%
% Manufactured solution U(t,x) = u0(sin(x^2 + omt) + eps)                 %
% dU/dt = om*u0*cos(x^2 + om*t)                                           %
% dU/dx = 2u0xcos(x^2 + om*t)                                             %
% d/dx(U^2/2) = 2u0^2*xcos(x^2 + om*t)(sin(x^2 + om*t) + eps)             %
% d^2U/dx^2 = 2u0(cos(x^2 + om*t) - 2x^2sin(x^2 + om*t))                  %
%-------------------------------------------------------------------------%

iter = 5;
const = 2.0;

for kk=1:iter
    
    dx = (xmax - xmin)/(npts - 1);
    xax = linspace(xmin,xmax,npts);
    
    % initial conditions
    source = zeros(1,npts);
    vx_new = u0*sin(xax.^2 + eps)*const;
    
    for ii=1
        
        % exact solution (steady)
        ex_sol = u0*sin(xax.^2) + eps;
        
        vx = vx_new;

        for jj=2:npts-1

            source(1,jj) = 2.0*u0^2*xax(jj)*cos(xax(jj)^2)*(sin(xax(jj)^2) + eps) - nu*(2.0*u0*(cos(xax(jj)^2) -...
                2.0*xax(jj)^2*sin(xax(jj)^2)));

            vx_new(1,jj) = vx(1,jj) - (dt/(2.0*dx))*((vx(1,jj+1)^2)/2 - (vx(1,jj-1)^2)/2) +...
                (dt/(dx^2))*nu*(vx(1,jj+1) - 2.0*vx(1,jj) + vx(1,jj-1)) + dt*source(1,jj);
        end

        % boundary conditions
        source(1,end) = c*cos(xmax + c*dt*ii) + (a + sin(xmax + ii*dt*c))*cos(xmax + c*ii*dt) + nu*sin(xmax + c*ii*dt);
        source(1,1) = c*cos(xmin + c*dt*ii) + (a + sin(xmin + ii*dt*c))*cos(xmin + c*ii*dt) + nu*sin(xmin + c*ii*dt);
%         vx_new(1,end) = vx(1,end) - (dt/(2.0*dx))*vx(1,end)*(vx(1,2) - vx(1,end-1)) +...
%                 (dt/(dx^2))*nu*(vx(1,2) - 2.0*vx(1,end) + vx(1,end-1)) + dt*source(1,end);
%         vx_new(1,1) = vx_new(1,end);

        vx_new(1,end) = a + sin(xmax + c*dt*ii);
        vx_new(1,1) = a + sin(xmin + c*dt*ii);
%        
        l_inf(1,kk) = norm(ex_sol - vx_new);
        l_sec(1,kk) = rms(ex_sol - vx_new);

        stop = 1;
        
        
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
