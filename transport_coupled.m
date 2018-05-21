%-------------------------------------------------------%
% solve coupled transport equations                     %
% continuity and momentum eqns                          %
% CONSERVATIVE FORMS                                    %
% (partial derivatives)                                 %
% dvz/dt + d(vz^2/2)/dz + (1/mn)(Te+Ti)dn/dz = 0        %
% dn/dt + d(n*vz)/dz = 0                                %
% staggered n and vz grids                              %
% momentum eqn central differenced                      %
% continuity eqn first order upwind (flux selecting)    %
% ghost points on density for mom source term           %
%       -- first order neumann, zero flux               %
% rlbarnett c3149416 140518                             %
%-------------------------------------------------------%

%------
% parameters %
%------
transport_large;

%%

figure(1)
set(gcf,'Position',[0   536   824   419])
plot(nxax,n_new,'DisplayName','time = 0s')
xlabel('Position (m)')
ylabel('Density m$^{-3}$')
xlim([min(nxax) max(nxax)])
drawnow
hold on

figure(2)
set(gcf,'Position',[859   536   824   419])
plot(vxax,vx_new,'DisplayName','time = 0s')
xlabel('Position (m)')
ylabel('Velocity ms$^{-1}$')
xlim([min(vxax) max(vxax)])
drawnow
hold on

count = 1;
timerVal = tic;

%-------------------------------------------------------------------------%
% Solve                                                                   %
%-------------------------------------------------------------------------%

for ii=1:nmax
    
    alpha = dt/(2.0*dx);
    
    vx = vx_new;
    n = n_new;
    n_interp = interp1(nxax,n,vxax,'spline');
    gradn = (n_interp(3:end) - n_interp(1:end-2))/(2.0*dx);

    for jj=2:npts-1
        
%---------------------------------------------------------------------------%
% central difference                                                        %
% currently unstable on staggered grid 070518                               %
%---------------------------------------------------------------------------%
%         nA(jj,jj) = 1.0 - alpha*(vx(1,jj+1) - vx(1,jj-1));
%         nA(jj,jj-1) = alpha*vx(1,jj);
%         nA(jj,jj+1) = -alpha*vx(1,jj);  
%---------------------------------------------------------------------------%
% first order upwind                                                        %
% stable on staggered grid                                                  %
%---------------------------------------------------------------------------%
        if ((vx(1,jj-1)+vx(1,jj))/2)>0
            nA(jj,jj) = 1.0 - alpha*vx(1,jj);
            nA(jj,jj-1) = alpha*vx(1,jj-1);
        elseif ((vx(1,jj-1)+vx(1,jj))/2)<=0
            nA(jj,jj) = 1.0 + alpha*vx(1,jj-1);
            nA(jj,jj+1) = -alpha*vx(1,jj);
            nA(2,2) = 1.0 + alpha*vx(1,1);
            nA(2,3) = -alpha*vx(1,2);
        end 
    end
    
    for jj=2:npts-2
        
        vxA(jj,jj) = (2.0*nu*dt)/(dx^2) - 1;
        vxA(jj,jj-1) = -(vx(1,jj-1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        vxA(jj,jj+1) = (vx(1,jj+1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        
        vx_source(jj,1) = -((Te + Ti)*e)/(m*n_interp(1,jj))*(gradn(1,jj-1)); 

    end
    
    nA = sparse(nA);
    vxA = sparse(vxA);
    
%     vx_source(1,1) = -((Te + Ti)*e)/(m*n_interp(1,1))*(gradn(1,1));
%     vx_source(end,1) = -((Te + Ti)*e)/(m*n_interp(1,end))*(gradn(1,end)); 
    
    n_new = nA*n';
    vx_new = dt*vx_source - vxA*vx';
    
    n_new = n_new';
    vx_new = vx_new';
    
    vx_new(1,1) = -cs;
    vx_new(1,end) = cs;
    % set first order neumann boundary conditions for the density ghost
    % points
    n_new(1,1) = n_new(1,2);
    n_new(1,end) = n_new(1,end-1);
    
    nan_check = isnan(vx_new);
    
    if (0.99*(dx^2)/(2.0*nu))<(0.99*dx/max(abs(vx_new)))
        dt = 0.99*(dx^2)/(2.0*nu);
    elseif (0.99*(dx^2)/(2.0*nu))>(0.99*dx/max(abs(vx_new)))
        dt = 0.99*dx/max(abs(vx_new));
    end
    
    if sum(nan_check) ~= 0
        fprintf('unstable, ii=%d\n',ii)
        break
    elseif dt*max(abs(vx_new))/dx >= 1.0 || dt*2*nu/dx^2 >= 1.0
        fprintf('CFL condition violated, ii=%d\n',ii)
        break
    elseif ii==count*round(nmax/10)
        fprintf('***--------------------***\n')
        fprintf('ii=%d, count=%d\n', [ii count])
        fprintf('dt=%ds\n', dt)
        fprintf('total time %d\n', toc(timerVal))
        if dt == 0.99*(dx^2)/(2.0*nu)
            fprintf('Diffusive CFL condition\n')
        elseif dt == 0.99*dx/max(abs(vx_new))
            fprintf('Convective CFL condition\n')
        end
        figure(1)
        plot(nxax,n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
%         ylim([-1.0e19 1.0e19]) 
        legend('show')
        drawnow
        hold on
%         pause(1)
        figure(2)
        plot(vxax,vx_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        legend('show')
        drawnow
        hold on
        figure(3)
        plot(vxax,vx_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        count = count + 1;
    end 
    
end

%%

figure(1)
legend('show','Location','northeast')
xlabel('Position (m)')
ylabel('Density m$^{-3}$')
hold off
drawnow

figure(2)
legend('show','Location','southeast')
xlabel('Position (m)')
ylabel('Velocity ms$^{-1}$')
hold off
drawnow

figure(3)
legend('show','Location','northwest')
xlabel('Position (m)')
ylabel('Velocity source ms$^{-1}$')
hold off





