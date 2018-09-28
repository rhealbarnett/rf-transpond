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
% m = mp;
m = 0.01;

%------
% parameters %
%------
Te = (1.0/e)*50.0;
Ti = (1.0/e)*50.0;
cs = 50.0;
nu = 0.2;


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

%------
% temporal domain %
%------
tmin = 0;
nmax = 1.0e7;

%%

%-------------------------------------------------------------------------%
% Set initial and boundary values for n and v                             %
% Initialise coefficient matrices                                         %
% calculate source term
%-------------------------------------------------------------------------%

u0 = 1.0;
n0 = 2.0;
om = 1.0e3;
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


%%

%-------------------------------------------------------------------------%
% Set tolerance; begin calculation                                        %
%-------------------------------------------------------------------------%

tol = 1.0e-14;
iter = 8;
l_inf = zeros(1,iter);
l_two = zeros(1,iter);
dx_arr = zeros(1,iter);


for kk=1:iter

    fprintf('counter=%d\n',kk)

    dx = (xmax - xmin)/(npts - 1);
    dx_arr(1,kk) = dx;
    xax = linspace(xmin,xmax,npts);
    mult = 1.0/dx;

    fprintf('dx=%d\n',dx)

     % initial conditions
    vx_new = u0*sin(xax.^2 + epsilon);
    n_new = n0*cos(xax.^2);
    rhs = gradient(n_new);
    vx_pos = sparse(npts,npts);
    vx_neg = sparse(npts,npts);
    vx_diff = sparse(npts,npts);
    vx_I = eye(npts,npts);
     
    for ii=1:nmax

        ex_sol = u0*(sin(xax.^2) + epsilon);

        vx = vx_new;

        for jj=2:npts-1
            if vx(1,jj)>0
                vx_pos(jj,jj) = - mult*vx(1,jj);
                vx_pos(jj,jj-1) = mult*vx(1,jj);
            elseif vx(1,jj)<0
                vx_neg(jj,jj) = mult*vx(1,jj);
                vx_neg(jj,jj+1) = -mult*vx(1,jj);
            end
            vx_diff(jj,jj) = - mult*((2.0*nu)/dx);
            vx_diff(jj,jj-1) = mult*(nu/dx);
            vx_diff(jj,jj+1) = mult*(nu/dx);
        end

        vxA = vx_I + vx_pos + vx_neg + vx_diff;
        vxA(1,1) = 1.0; vxA(end,end) = 1.0;

        source_dx = 2.0*u0*xax.*cos(xax.^2);
        source_dxx = 2.0*u0*(cos(xax.^2) - 2.0*xax.^2.*sin(xax.^2));
        source_nx = -2.0*n0*xax.*sin(xax.^2);

        source = vx.*source_dx + nu*source_dxx - source_nx;

        vx_new = vxA\(source' + rhs');

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

    npts = 2.0*npts-1;
     
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

























