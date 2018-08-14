%-------------------------------------------------------%
% solve coupled transport equations                     %
% continuity and momentum eqns                          %
% CONSERVATIVE FORM OF MASS CONT.                       %
% (partial derivatives)                                 %
% dvz/dt + d(vz^2/2)/dz + (1/mn)(Te+Ti)dn/dz = 0        %
% dn/dt + d(n*vz)/dz = 0                                %
% colocated n and vz grids                              %
% momentum eqn central differenced                      %
% continuity eqn first order upwind (flux selecting)    %
% ghost points on density for mom source term           %
%       -- first order neumann, zero flux               %
% rlbarnett c3149416 140518                             %
%-------------------------------------------------------%

%------
% parameters %
%------
transport_params_colocated;

%%

figure(1)
set(gcf,'Position',[563 925 560 420])
plot(nxax(2:npts-1),n_new(2:npts-1),'DisplayName','time = 0s')
hold on

figure(2)
set(gcf,'Position',[7 925 560 420])
plot(vxax,vx_new/cs,'DisplayName','time = 0s')
hold on

figure(3)
set(gcf,'Position',[3 476 560 420])
plot(vxax,vx_source*dt,'DisplayName','time = 0s')
hold on

count = 1;
timerVal = tic;

vx_mat = zeros(nmax,npts-1);
n_mat = zeros(nmax,npts-2);
pressure_mat = zeros(nmax,npts-2);

%%
%-------------------------------------------------------------------------%
% Solve                                                                   %
%-------------------------------------------------------------------------%

for ii=1:nmax
    
    alpha = dt/dx;
    
    vx = vx_new;
    n = n_new;

    for jj=2:npts-1
%---------------------------------------------------------------------------%
% first order upwind                                                        %
% stable on staggered grid                                                  %
%---------------------------------------------------------------------------%
        if vx(1,jj)>0
            nA(jj,jj) = 1.0 - alpha*vx(1,jj);
            nA(jj,jj-1) = alpha*vx(1,jj-1);
            nA(npts,npts) = 1.0 - alpha*vx(1,npts);
            nA(npts,npts-1) = alpha*vx(1,npts-1);
        elseif vx(1,jj)<0
            nA(jj,jj) = 1.0 + alpha*vx(1,jj);
            nA(jj,jj+1) = -alpha*vx(1,jj+1);
            nA(1,1) = 1.0 + alpha*vx(1,1);
            nA(1,2) = -alpha*vx(1,2);
        end    
    end
    
    for jj=2:npts-1
        
        if vx(1,jj)>0
            vxA(jj,jj) = 1 - alpha*vx(1,jj);
            vxA(jj,jj-1) = alpha*vx(1,jj);
            vxA(end,end) = 1 - alpha*vx(1,end);
            vxA(end,end-1) = alpha*vx(1,end);
            pond_source(jj,1) = (1.0/m)*pond_const*((Efield(1,jj) - Efield(1,jj-1))/dx);
            pond_source(npts,1) = (1.0/m)*pond_const*((Efield(1,npts) - Efield(1,npts-1))/dx);
            vx_source(jj,1) = -((Te + Ti)*e/(m*n(1,jj)))*((n(1,jj) - n(1,jj-1))/dx) -...
                pond_source(jj,1);
            vx_source(npts,1) = -((Te + Ti)*e/(m*n(1,npts)))*((n(1,npts) - n(1,npts-1))/dx) -...
                pond_source(npts,1);
        elseif vx(1,jj)<0
            vxA(jj,jj) = 1 + alpha*vx(1,jj);
            vxA(jj,jj+1) = -alpha*vx(1,jj);
            vxA(1,1) = 1 + alpha*vx(1,1);
            vxA(1,2) = -alpha*vx(1,1);
            pond_source(jj,1) = (1.0/m)*pond_const*((Efield(1,jj+1) - Efield(1,jj))/dx);
            pond_source(1,1) = (1.0/m)*pond_const*((Efield(1,2) - Efield(1,1))/dx);
            vx_source(jj,1) = -((Te + Ti)*e/(m*n(1,jj)))*((n(1,jj+1) - n(1,jj))/dx) -...
                pond_source(jj,1);
            vx_source(1,1) = -((Te + Ti)*e/(m*n(1,1)))*((n(1,2) - n(1,1))/dx) -...
                pond_source(1,1);
        end

    end

    nA = sparse(nA);
    vxA = sparse(vxA);
   
    vx_new = dt*vx_source + vxA*vx';
    n_new = dt*n_source + nA*n';
    
    vx_new = vx_new';
    n_new = n_new';
    
%     if vx(1,jj)>0
%         vx_new(1,1) = vx_new(1,2);
%     elseif vx(1,jj)<0
%         vx_new(1,npts-1) = vx_new(1,npts-2);
%     end

    vx_new(1,1) = vx_new(1,2);
%     n_new(1,1) = (4.0*n_new(1,2) - n_new(1,3))/3.0;
    n_new(1,1) = n_new(1,2);
    
    nan_check = isnan(vx_new);
    
    if (cfl_fact*(dx^2)/(2.0*nu))<(cfl_fact*dx/max(abs(vx_new)))
        dt = cfl_fact*(dx^2)/(2.0*nu);
    elseif (cfl_fact*(dx^2)/(2.0*nu))>(cfl_fact*dx/max(abs(vx_new)))
        dt = cfl_fact*dx/max(abs(vx_new));
    end
    
    if sum(nan_check) ~= 0
        fprintf('unstable, ii=%d\n',ii)
        break
    elseif dt*max(abs(vx_new))/dx >= 1.0 || dt*2*nu/dx^2 >= 1.0
        fprintf('CFL condition violated, ii=%d\n',ii)
        break
    elseif ii==count*round(nmax/5)
        fprintf('***--------------------***\n')
        fprintf('ii=%d, count=%d\n', [ii count])
        fprintf('dt=%ds\n', dt)
        fprintf('total time=%ds\n', dt*ii)
        fprintf('simulation time %d\n', toc(timerVal))
        if dt == 0.8*(dx^2)/(2.0*nu)
            fprintf('Diffusive CFL condition\n')
        elseif dt == 0.8*dx/max(abs(vx_new))
            fprintf('Convective CFL condition\n')
        end
        figure(1)
        set(gcf,'Position',[563 925 560 420])
        plot(nxax,n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax) max(nxax)])
        hold on
        figure(2)
        set(gcf,'Position',[7 925 560 420])
        plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(3)
        set(gcf,'Position',[3 476 560 420])
        plot(vxax,vx_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(6)
        set(gcf,'Position',[3 33 560 420])
        plot(vxax,pond_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax) max(nxax)])
        hold on
        count = count + 1;
    end 
    
end

%%

figure(1)
set(gcf,'Position',[563 925 560 420])
plot(nxax,n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density m^{-3}','Fontsize',16)
legend('show','Location','south')
hold off

figure(2)
set(gcf,'Position',[7 925 560 420])
plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Mach number','Fontsize',16)
legend('show','Location','southeast')
hold off

figure(3)
set(gcf,'Position',[3 476 560 420])
plot(vxax,vx_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Velocity source ms^{-1}','Fontsize',16)
legend('show','Location','northwest')
hold off

figure(6)
set(gcf,'Position',[3 33 560 420])
plot(vxax,pond_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Ponderomotive source ms^{-1}','Fontsize',16)
legend('show','Location','south')
hold off

