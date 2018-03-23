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
dt = 1.e-5;

% max time to run
tmax = 1.0e4*dt;

% max number of iterations -- keep constant
nmax = tmax/dt;

%------
% spatial domain %
%------
xmin = 0;
xmax = 1;

% vary spatial step size
dx = logspace(-1,-5,5);

%%
%-------------------------------------------------------------------------%
% Manufactured solution U(t,x) = sin(x-t)                                 %
% dU/dt = -cos(x-t); d/dx(U^2/2) = -(1/2)sin(2(t-x))                      %
%-------------------------------------------------------------------------%

for kk=1
    kk
    for ii=1:nmax
        
        % number of spatial points will depend on dx
        npts = round((xmax-xmin)/dx(kk));
        xax = linspace(xmin,xmax,npts);
        
        % initial conditions
        vx = sin(xax);
        
        % boundary conditions (Dirichlet, time dependent)
        vx(1,1) = sin(-dt*ii);
        vx(1,end) = sin(xmax - dt*ii);
        
        % initialise source, error arrays
        source = zeros(1,npts);
        l_inf = zeros(1,npts);
        l_sec = zeros(1,npts);
        
        % source boundary conditions
        source(1,1) = - cos(ii*dt) - (1.0/2.0)*sin(2.0*(ii*dt));
        source(1,end) = - cos(ii*dt - xmax) - (1.0/2.0)*sin(2.0*(ii*dt - xmax));
        
        for jj=2:npts-1
            
            % discretised solution
            vx(1,jj) = (1./2.)*(vx(1,jj-1) + vx(1,jj+1)) - (dt/(2.0*dx(kk)))*((vx(1,jj+1)^2)/2 - (vx(1,jj-1)^2/2));
            
            % 'exact' solution
            source(1,jj) = - cos(ii*dt - dx(kk)*jj) - (1.0/2.0)*sin(2.0*(ii*dt - dx(kk)*jj)); 

        end

        l_inf(1,ii) = abs(max(vx - source));
        l_sec(1,ii) = rms(vx - source);
        
    end
end


%%
%------ 
% plot solution %
%------

figure(1)
plot(xax, vx)
hold on
plot(xax, source)
legend('vx', 'source')
hold off
