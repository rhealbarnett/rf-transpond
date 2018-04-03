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

% max number of iterations -- keep constant
nmax = 1.0e7;

% temporal step
% dt = 1.0/nmax;
dt = 0.0;

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
epsilon = 0.001;
nu = 0.7;

%%
%-------------------------------------------------------------------------%
% Manufactured solution U(t,x) = u0(sin(x^2 + omt) + eps)                 %
% dU/dt = om*u0*cos(x^2 + om*t)                                           %
% dU/dx = 2u0xcos(x^2 + om*t)                                             %
% d/dx(U^2/2) = 2u0^2*xcos(x^2 + om*t)(sin(x^2 + om*t) + eps)             %
% d^2U/dx^2 = 2u0(cos(x^2 + om*t) - 2x^2sin(x^2 + om*t))                  %
%-------------------------------------------------------------------------%

iter = 8;
const = 2.0;
tol = 1.0e-14;
l_inf = zeros(1,iter);
l_two = zeros(1,iter);
dx_arr = zeros(1,iter);

om = 0.0;


for kk=1:iter
    
    fprintf('counter=%d\n',kk)
    
    dx = (xmax - xmin)/(npts - 1);
    dx_arr(1,kk) = dx;
    xax = linspace(xmin,xmax,npts);
    
    fprintf('dx=%d\n',dx)

    % initial conditions
    source = zeros(npts,1);
    vx_new = u0*sin(xax.^2 + epsilon)*const;
    coeff_mat = zeros(npts,npts);
    
    for ii=1:nmax
        
        % exact solution (steady)
        ex_sol = u0*sin(xax.^2) + epsilon;
        
        % boundary conditions
        coeff_mat(1,1) = 1.0;
        coeff_mat(end,end) = 1.0;
        
        source(1,1) = u0*sin(xmin.^2) + epsilon;
        source(end,1) = u0*sin(xmax.^2) + epsilon;
        
        vx = vx_new;

        for jj=2:npts-1
            
            coeff_mat(jj,jj) = (2.0*nu)/(dx^2);
            coeff_mat(jj,jj-1) = -vx(1,jj)/(2.0*dx) - nu/(dx^2);
            coeff_mat(jj,jj+1) = vx(1,jj)/(2.0*dx) - nu/(dx^2);
            
            coeff_mat = sparse(coeff_mat);

            source(jj,1) = 2.0*u0^2*xax(1,jj)*cos(xax(1,jj)^2)*(sin(xax(1,jj)^2) +...
                epsilon) - nu*2.0*u0*(cos(xax(1,jj)^2) - 2.0*xax(1,jj)^2*sin(xax(1,jj)^2));

%             vx_new(1,jj) = (1.0/2.0)*(vx(1,jj+1) - vx(1,jj-1)) - (dx/(4.0*nu))*(vx(1,jj+1)^2/2.0 - vx(1,jj-1)^2/2.0) +...
%                 (dx^2/(2.0*nu))*source(1,jj);
            
        end
        
        
        vx_new = coeff_mat\source;
        
        vx_new = vx_new';
        
        l_inf(1,kk) = norm(ex_sol - vx_new, Inf);
        l_two(1,kk) = rms(ex_sol - vx_new);
        
        if rms(vx_new - vx)<=tol
            fprintf('tolerance reached, ii=%d\n',ii)
            break
        else
            continue
        end
        
    end
    
    npts = 2.0*npts - 1;
    fprintf('coefficient matrix determinant=%d\n\n',det(coeff_mat))
    
end

ratio_inf = l_inf(1:iter-1)./l_inf(2:iter);
ratio_two = l_two(1:iter-1)./l_two(2:iter);

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
loglog(dx_arr,l_two, '-*r') 
hold on
loglog(dx_arr,l_inf, '-ob')
loglog(dx_arr, dx_arr.^2, '--k')
% xlim([5e-4 1e-1])
xlabel('dx','Fontsize',16)
ylabel('Error','Fontsize',16)
legend({'L$_2$', 'L$_{\inf}$', 'dx$^2$'},'Fontsize',16,'Location','northwest')
hold off
