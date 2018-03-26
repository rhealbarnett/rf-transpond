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

% spatial step and axis
dx = logspace(-5,-1,5);

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

l_sec = zeros(1,length(dx));
l_inf = zeros(1,length(dx));

for kk=1:length(dx)
    
    npts = round(1 + (xmax - xmin)/dx(kk));
    xax = linspace(xmin,xmax,npts);
    
    % initial conditions
    source = zeros(1,npts);
    vx = zeros(1,npts);
    vx_new = a + sin(xax) + source;
    
    for ii=1

        % exact solution
        ex_sol = a + sin(xax + c*dt*ii);

        vx = vx_new;

        for jj=2:npts-1

            source(1,jj) = c*cos(c*ii*dt + xax(jj)) + (a + sin(xax(jj) + ii*dt*c))*cos(xax(jj) + c*ii*dt) + nu*sin(xax(jj) + c*ii*dt); 

            vx_new(1,jj) = vx(1,jj) - ((dt)/(2.0*dx(kk)))*((vx(1,jj+1)^2)/2 - (vx(1,jj-1)^2)/2) +...
                ((nu*dt)/(dx(kk)^2))*(vx(1,jj+1) - 2.0*vx(1,jj) + vx(1,jj-1)) + dt*source(1,jj);
        end

        % boundary conditions
        vx_new(1,end) = a + sin(xmax + c*dt*ii);
        vx_new(1,1) = a + sin(c*dt*ii);
        
        l_inf(1,kk) = abs(max(vx_new - ex_sol));
        l_sec(1,kk) = rms(vx_new - ex_sol);

    end
    
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
