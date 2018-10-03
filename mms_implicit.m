%--------------------------------------------------------------------------------------------------------------%
% MMS for implicit solve on momentum equation                                                                            %
% rlbarnett c3149416 280918                                                                                    %
%--------------------------------------------------------------------------------------------------------------%
% u(x,t) = u0*(sin(x^2 + om*t) + epsilon)
% dtu(x,t) = u0*om*cos(x^2 + om*t)
% dxu(x,t) = 2*u0*x*cos(x^2 + om*t)
% u(x,t)*dxu(x,t) = 2*u0^2*x*(sin(x^2 + om*t) + epsilon)*cos(x^2 + om*t)
% dxxu(x,t) = 2*u0*(cos(x^2 + om*t) - 2*x^2*sin(x^2 + om*t))
% n(x) = n0*cos(x^2)
% dxn(x) = -2*n0*x*sin(x^2)
%--------------------------------------------------------------------------------------------------------------%

%------
% constants %
%------
constants;
e = const.e;
% m = mp;
m = 0.01;

%------
% parameters %
%------
Te = (1.0/e)*50.0;
Ti = (1.0/e)*50.0;
cs = 50.0;

%------
% spatial domain %
%------
xmin = -0.03;
xmax = 0.1;

% include two additional gridpoints for the density ghost points
% velocity grid will then be defined as having npts-1 (xax(1:npts-1)) --
% density solution space will be defined as having npts-2 (xax(2:npts-1))
npts = 11;
xax = linspace(xmin,xmax,npts);
dx = (xmax - xmin)/(npts - 1);

%------
% temporal domain %
%------
tmin = 0;

% dt = 8.0e-6;
nmax = 1.0e7;
dt = 1.0;


%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
% calculate source term
%-------------------------------------------------------------------------%

u0 = 1.0;
n0 = 2.0;
nu = 0.2;
% om = 1.0e3;
om = 0.0;
epsilon = 0.001;

vx_new = u0*(sin(xax.^2) + epsilon);
n_new = n0*cos(xax.^2);
rhs = gradient(n_new);

beta = -((Te + Ti)*e./(m*n_new));
rhs = beta.*gradient(n_new);

vx_pos = sparse(npts,npts);
vx_neg = sparse(npts,npts);
vx_diff = sparse(npts,npts);
vx_I = eye(npts,npts);
source_dt = zeros(1,npts);
source_dx = zeros(1,npts);
source_dxx = zeros(1,npts);


%%

%-------------------------------------------------------------------------%
% Set tolerance; begin calculation                                        %
%-------------------------------------------------------------------------%

tol = 1.0e-14;
iter = 8;
r = 2.0;
count = 1;
l_inf = zeros(1,iter);
l_two = zeros(1,iter);
dx_arr = zeros(1,iter);
dt_arr = zeros(1,iter);

for kk=1:iter

    fprintf('counter=%d\n',kk)

    dx = (xmax - xmin)/(npts - 1);
    dx_arr(1,kk) = dx;
%     tmax = 4.0e-4;
%     nmax = round(tmax/dt);
%     dt_arr(1,kk) = dt;
    xax = linspace(xmin,xmax,npts);
    fprintf('dx=%d\n',dx)

     % initial conditions
    vx_new = u0*(sin(xax.^2) + epsilon);
    n_new = n0*cos(xax.^2);
    rhs = gradient(n_new);
    vx_pos = sparse(npts,npts);
    vx_neg = sparse(npts,npts);
    vx_diff = sparse(npts,npts);
    vx_I = eye(npts,npts);
    source_dt = zeros(1,npts);
    source_dx = zeros(1,npts);
    source_dxx = zeros(1,npts);
     
    for ii=1:nmax

        ex_sol = u0*(sin(xax.^2 + om*dt*ii) + epsilon);

        vx = vx_new;
       
        for jj=2:npts-1
            if vx(1,jj)>0
                vx_pos(jj,jj) = - vx(1,jj)/dx;
                vx_pos(jj,jj-1) = vx(1,jj)/dx;
            elseif vx(1,jj)<0
                vx_neg(jj,jj) = vx(1,jj)/dx;
                vx_neg(jj,jj+1) = -vx(1,jj)/dx;
            end
            vx_diff(jj,jj) = - (2.0*nu)/(dx^2);
            vx_diff(jj,jj-1) = (nu/(dx^2));
            vx_diff(jj,jj+1) = (nu/(dx^2));
        end

        vxA = vx_pos + vx_neg + vx_diff;
        Avx = vx_I - dt*vxA;
%         Avx = vx_I + dt*vxA;
        Avx(1,1) = 1.0; Avx(end,end) = 1.0;
        vx(1,1) = u0*(sin(xmin^2 + om*dt*ii) + epsilon);
        vx(1,end) = u0*(sin(xmax^2 + om*dt*ii) + epsilon);
    
        source_dt = u0*om*cos(xax.^2 + om*dt*ii);
        source_dx = 2.0*u0*xax.*cos(xax.^2 + om*dt*ii)*u0.*(sin(xax.^2 + om*dt*ii) + epsilon);
        source_dxx = 2.0*u0*(cos(xax.^2 + om*dt*ii) - 2.0*xax.^2.*sin(xax.^2 + om*dt*ii));
        source_nx = -2.0*n0*xax.*sin(xax.^2);

        source = source_dt + source_dx - nu*source_dxx;
        
        vx_new = Avx\(source');    
%         vx_new = Avx\(vx' + dt*source');
%         vx_new = Avx*vx' + dt*source';

        vx_new = vx_new';

        l_inf(1,kk) = norm(ex_sol - vx_new, Inf);
        l_two(1,kk) = rms(ex_sol - vx_new);

        if rms(vx_new - vx)<=tol
            fprintf('tolerance reached, ii=%d\n',ii)
            break
%         elseif ii==count*100
%             fprintf('ii=%d\n',ii)
%             fprintf('rms=%d\n',rms(vx_new - vx))
%             count = count + 1;
%             continue
        else 
            continue
        end
        
        
    end

    npts = r*npts-1;
%     dt = dt/r;
     
end

ratio_inf = l_inf(1:iter-1)./l_inf(2:iter);
ratio_two = l_two(1:iter-1)./l_two(2:iter);
oo_inf = log(ratio_inf)/log(2);
oo_two = log(ratio_two)/log(2);

%%
%------ 
% plot solution %
%------

figure(1)
plot(xax, vx, '*k')
hold on
plot(xax, ex_sol)
% xlabel('Location')
% ylabel('Amplitude')
legend('vx', 'exact solution')
hold off

figure(2)
loglog(dx_arr,l_two, '-*r') 
hold on
loglog(dx_arr,l_inf, '-ob')
loglog(dx_arr, dx_arr, '--k')
% xlim([5e-4 1e-1])
% xlabel('dx')
% ylabel('Error')
legend({'L_2', 'L_{\infty}', 'dx'},'Location','northwest')
hold off

























