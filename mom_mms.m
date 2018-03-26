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
xmax = 1.0;

% number of spatial points
npts = 61;

% spatial step and axis
dx = (xmax - xmin)/(npts-1);
xax = linspace(xmin,xmax,npts);

% coefficients for source term s = a + sin(x + ct)
a = 1.0;
c = 1.0;

%%
%-------------------------------------------------------------------------%
% Manufactured solution U(t,x) = sin(x-t)                                 %
% dU/dt = -cos(x-t); d/dx(U^2/2) = -(1/2)sin(2(t-x))                      %
%-------------------------------------------------------------------------%

% diffusion term multiplier
nu = 0.0;

% initial conditions
source = zeros(1,npts);%c*cos(xax) + (a + sin(xax)).*cos(xax) + nu*sin(xax);
vx = zeros(1,npts);
vx_new = a + sin(xax) + source;
% vx_new = zeros(1,npts);

l_sec = zeros(1,nmax);
l_inf = zeros(1,nmax);

for ii=1
    
    % exact solution
    ex_sol = a + sin(xax + c*dt*ii);
    
    vx = vx_new;

    for jj=2:npts-1

        source(1,jj) = c*cos(c*ii*dt + xax(jj)) + (a + sin(xax(jj) + ii*dt*c))*cos(xax(jj) + c*ii*dt) + nu*sin(xax(jj) + c*ii*dt); 

        vx_new(1,jj) = vx(1,jj) - ((dt)/(2.0*dx))*((vx(1,jj+1)^2)/2 - (vx(1,jj-1)^2)/2) +...
            ((nu*dt)/(dx^2))*(vx(1,jj+1) - 2.0*vx(1,jj) + vx(1,jj-1)) + dt*source(1,jj);
    end

    l_inf(1,ii) = abs(max(vx_new - ex_sol));
    l_sec(1,ii) = rms(vx_new - ex_sol);
    
    % source boundary conditions
    source(1,1) = c*cos(ii*dt*c) + (a + sin(c*dt*ii))*cos(c*dt*ii) + nu*sin(c*ii*dt);
    source(1,end) = c*cos(ii*dt*c + xmax) + (a + sin(xmax + c*dt*ii))*cos(xmax + c*dt*ii) + nu*sin(xmax + c*ii*dt);

%     vx_new(1,end) = vx(1,end) - ((dt)/(2.0*dx))*((vx(1,2)^2)/2 - (vx(1,end-1)^2)/2) +...
%             ((nu*dt)/(dx^2))*(vx(1,2) - 2.0*vx(1,end) + vx(1,end-1)) + dt*source(1,end);
    vx_new(1,end) = a + sin(xmax + c*dt*ii);
%     vx_new(1,1) = vx_new(1,end);
    vx_new(1,1) = a + sin(c*dt*ii);

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
xlabel('Location')
ylabel('Amplitude')
legend('vx', 'exact solution')
hold off
