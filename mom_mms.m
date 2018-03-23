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

% choose dt based on smallest dx (assume CFL condition is dt/dx <= 1)
dt = 5.e-4*0.99;

% max time to run
tmax = 1.0e4*dt;

% max number of iterations -- keep constant
nmax = tmax/dt;

%------
% spatial domain %
%------
xmin = 0;
xmax = 2.0*pi;

% vary spatial step size
dx = logspace(-1,-5,5);

% coefficients for source term s = a + sin(x + ct)
a = 1.0;
c = 1.0;

%%
%-------------------------------------------------------------------------%
% Manufactured solution U(t,x) = sin(x-t)                                 %
% dU/dt = -cos(x-t); d/dx(U^2/2) = -(1/2)sin(2(t-x))                      %
%-------------------------------------------------------------------------%

nu = 0.05;
dx = 0.1047;
% for kk=length(dx)
%     kk
    
% number of spatial points will depend on dx
npts = round((xmax-xmin)/dx);
xax = linspace(xmin,xmax,npts);

vx_new = zeros(1,npts);

% initial conditions
source = 0.0;%c*cos(xax) + (a + sin(xax)).*cos(xax) + nu*sin(xax);
vx = sin(xax) + source;

for ii=1:nmax


    % source boundary conditions
%         source(1,1) = - cos(ii*dt) + sin(-ii*dt)*cos(-ii*dt) + nu*sin(-ii*dt);
    source(1,end) = c*cos(ii*dt*c + xmax) + (a + sin(xmax + c*dt*ii))*cos(xmax + c*dt*ii) + nu*sin(xmax + c*ii*dt);

    % boundary conditions (Dirichlet, time dependent)
%         vx(1,1) = sin(-dt*ii*c) + dt*source(1,1);
%     vx(1,end) = vx(1,end) - ((dt*vx(1,end))/(2.0*dx))*(vx(1,2) - vx(1,end-1)) +...
%             ((nu*dt)/(dx^2))*(vx(1,1) - 2.0*vx(1,end) + vx(1,end-1)) + dt*source(1,end);
% 
%     vx(1,1) = vx(1,end);

    % initialise source, error arrays
    l_inf = zeros(1,npts);
    l_sec = zeros(1,npts);

    % exact solution
    ex_sol = a + sin(xax + c*dt*ii);

    % calculation for second last spatial point
%         source(1,end-1) = - cos(ii*dt - xax(end-1)) + sin(xax(end-1) - ii*dt)*cos(xax(end-1) - ii*dt) + nu*sin(xax(end-1) - ii*dt); 

    % backward difference for second last spatial point (not included
    % in loop)
%         vx(1,end-1) = vx(1,end-1) - ((dt*vx(1,end-1))/(dx(kk)))*(vx(1,end-1) - vx(1,end-2)) +...
%             ((nu*dt)/(dx(kk)^2))*(vx(1,end-1) - 2.0*vx(1,end-2) + vx(1,end-3)) + dt*source(1,end-1);

    for jj=2:npts-1

        source(1,jj) = c*cos(c*ii*dt + xax(jj)) + (a + sin(xax(jj) + ii*dt*c))*cos(xax(jj) + c*ii*dt) + nu*sin(xax(jj) + c*ii*dt); 

        % discretised solution
        vx_new(1,jj) = vx(1,jj) - ((dt*vx(1,jj))/(2.0*dx))*(vx(1,jj+1) - vx(1,jj-1)) +...
            ((nu*dt)/(dx^2))*(vx(1,jj+1) - 2.0*vx(1,jj) + vx(1,jj-1)) + dt*source(1,jj);
    end

    l_inf(1,ii) = abs(max(vx - ex_sol));
    l_sec(1,ii) = rms(vx - ex_sol);

    vx_new(1,end) = vx(1,end) - ((dt*vx(1,end))/(2.0*dx))*(vx(1,2) - vx(1,end-1)) +...
            ((nu*dt)/(dx^2))*(vx(1,2) - 2.0*vx(1,end) + vx(1,end-1)) + dt*source(1,end);
    vx_new(1,1) = vx_new(1,end);    
    vx = vx_new;
end
% end


%%
%------ 
% plot solution %
%------

figure(1)
plot(xax, vx, '*k')
hold on
plot(xax, ex_sol)
legend('vx', 'exact solution')
hold off
