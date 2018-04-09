%----------------------------------------%
% Method of manufactured solutions       %
% mom eqn                                %
% dvx/dt + dvx/dx = 0                    %
% steady state solution                  %
% rlbarnett c3149416 300318              %
%----------------------------------------%

%------
% temporal domain %
%------
tmin = 0;

% temporal step
dt = 8.0e-6;
% dt = 1.0;

nmax = 1.0e7;

%------
% spatial domain %
%------
xmin = -0.1;
xmax = 1.7;

% initial grid size
npts = 1281;
dx = (xmax - xmin)/(npts - 1);
xax = linspace(xmin,xmax,npts);
% dx_arr(1,kk) = dx;

% coefficients for source term u0(sin(x^2 + omt) + eps) and viscosity
u0 = 1.0;
epsilon = 0.001;
nu = 0.7;
om = 1.0e3;

%%
%-------------------------------------------------------------------------%
% Manufactured solution U(t,x) = u0(sin(x^2 + omt) + eps)                 %
% dU/dt = om*u0*cos(x^2 + om*t)                                           %
% dU/dx = 2u0xcos(x^2 + om*t)                                             %
% d/dx(U^2/2) = 2u0^2*xcos(x^2 + om*t)(sin(x^2 + om*t) + eps)             %
% d^2U/dx^2 = 2u0(cos(x^2 + om*t) - 2x^2sin(x^2 + om*t))                  %
%-------------------------------------------------------------------------%

iter = 8;
const = 1.0;
tol = 1.0e-14;
l_inf = zeros(1,iter);
l_two = zeros(1,iter);
dx_arr = zeros(1,iter);
dt_arr = zeros(1,iter);


for kk=1:iter
    
    fprintf('counter=%d\n',kk)
    
    fprintf('dt=%d\n',dt)
  
    tmax = 4.0e-4;
    nmax = tmax/dt;
    dt_arr(1,kk) = dt;
    
    fprintf('nmax=%d\n',nmax)
    fprintf('ratio=%d\n',(dt/dx))

    % initial conditions
    source = zeros(npts,1);
    vx_new = u0*(sin(xax.^2) + epsilon)*const;
    coeff_mat = zeros(npts,npts);
    
    for ii=1:nmax
        
        % exact solution (steady)
        ex_sol = u0*(sin(xax.^2 + dt*ii*om) + epsilon);
        
        if ii==1
            figure(3)
            plot(xax, ex_sol)
            hold on
        end
        
        % boundary conditions
        coeff_mat(1,1) = 1.0;
        coeff_mat(end,end) = 1.0;
        
%         source(1,1) = u0*(sin(xmin.^2 + om*dt*ii) + epsilon);
%         source(end,1) = u0*(sin(xmax.^2 + om*dt*ii) + epsilon);
        
        vx = vx_new;

        for jj=2:npts-1
            
            coeff_mat(jj,jj) = (2.0*nu*dt)/(dx^2) - 1;
            coeff_mat(jj,jj-1) = -(vx(1,jj-1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
            coeff_mat(jj,jj+1) = (vx(1,jj+1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
            
            coeff_mat = sparse(coeff_mat);

            source(jj,1) = dt*(om*u0*cos(xax(1,jj)^2 + om*dt*ii) + 2.0*u0^2*xax(1,jj)*cos(xax(1,jj)^2 + om*dt*ii)*(sin(xax(1,jj)^2 + dt*om*ii) +...
                epsilon) - nu*2.0*u0*(cos(xax(1,jj)^2 + om*dt*ii) - 2.0*xax(1,jj)^2*sin(xax(1,jj)^2 + om*dt*ii)));
            
        end
        
        
        vx_new = source - coeff_mat*vx';
        
        vx_new(1,1) = u0*(sin(xmin.^2 + om*dt*ii) + epsilon);
        vx_new(end,1) = u0*(sin(xmax.^2 + om*dt*ii) + epsilon);
        
        vx_new = vx_new';
        
        l_inf(1,kk) = norm(ex_sol - vx_new, Inf);
        l_two(1,kk) = rms(ex_sol - vx_new);
        
%         if rms(vx_new - vx)<=tol
%             fprintf('tolerance reached, ii=%d\n',ii)
%             break
%         else
%             continue
%         end
        
    end
    
    dt = dt/2.0;
%     npts = 2*npts - 1;
    fprintf('coefficient matrix determinant=%d\n\n',det(coeff_mat))
    
end

ratio_inf = l_inf(1:iter-1)./l_inf(2:iter);
ratio_two = l_two(1:iter-1)./l_two(2:iter);

figure(3)
plot(xax, ex_sol, '*r')
legend('1 time step','100 time steps')
hold off

%%
%------ 
% plot solution %
%------

figure(1)
plot(xax, vx_new, '*k')
hold on
plot(xax, ex_sol)
xlabel('Location')
ylabel('Amplitude')
legend('vx', 'exact solution')
hold off

figure(2)
loglog(dt_arr, l_two, '-*r') 
hold on
loglog(dt_arr, l_inf, '-ob')
loglog(dt_arr, dt_arr, '--k')
% xlim([5e-4 1e-1])
xlabel('dt','Fontsize',16)
ylabel('Error','Fontsize',16)
legend({'L$_2$', 'L$_{\inf}$', 'dt'},'Fontsize',16,'Location','northwest')
hold off
