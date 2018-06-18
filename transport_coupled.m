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
semilogy(nxax,n_new,'DisplayName','time = 0s')
hold on

figure(2)
plot(vxax,vx_new/cs,'DisplayName','time = 0s')
hold on

figure(3)
plot(vxax,vx_source*dt,'DisplayName','time = 0s')
hold on

figure(4)
semilogy(nxax,n_source*dt,'DisplayName','time = 0s')
hold on

count = 1;
timerVal = tic;

%%
%-------------------------------------------------------------------------%
% Solve                                                                   %
%-------------------------------------------------------------------------%

for ii=1:nmax
    
    alpha = dt/(2.0*dx);
    
    vx = vx_new;
    n = n_new;
    n_fit = polyfit(nxax,n,2);
    n_interp = polyval(n_fit,vxax);
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
%             nA(2,2) = 1.0 + alpha*vx(1,1);
%             nA(2,3) = -alpha*vx(1,2);
        end 
        
        n_source(jj,1) = n(1,jj)*n_neut(jj,1)*rate_coeff;
        
    end
    
    n_source(1,1) = n_source(2,1);
    n_source(end,1) = n_source(end-1,1);
    
    for jj=2:npts-2
        
        vxA(jj,jj) = (2.0*nu*dt)/(dx^2) - 1;
        vxA(jj,jj-1) = -(vx(1,jj-1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        vxA(jj,jj+1) = (vx(1,jj+1)*dt)/(4.0*dx) - (nu*dt)/(dx^2);
        
        pond_source(jj,1) = 0.0;%(1.0/m)*(pond_pot(1,jj+1) - pond_pot(1,jj-1))/(2.0*dx);
        vx_source(jj,1) = -((Te + Ti)*e)/(m*n_interp(1,jj))*(gradn(1,jj-1)) -...
            pond_source(jj,1); 

    end
    
    nA = sparse(nA);
    vxA = sparse(vxA);
    
    n_new = dt*n_source + nA*n';
    vx_new = dt*vx_source - vxA*vx';
    
    n_new = n_new';
    vx_new = vx_new';
     
    vx_new(1,1) = -cs;
    vx_new(1,end) = cs;
    % set first order neumann boundary conditions for the density ghost
    % points
    n_new(1,1) = n_new(1,2);
    n_new(1,end) = n_new(1,end-1);
    
    l_inf_vx(1,ii) = norm(vx - vx_new);
    l_two_vx(1,ii) = rms(vx - vx_new);
    l_inf_n(1,ii) = norm(n - n_new)/norm(n);
    l_two_n(1,ii) = rms(n - n_new);
    
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
    elseif ii==count*round(nmax/6)
        fprintf('***--------------------***\n')
        fprintf('ii=%d, count=%d\n', [ii count])
        fprintf('dt=%ds\n', dt)
        fprintf('total time=%ds\n', dt*ii)
        fprintf('simulation time %d\n', toc(timerVal))
        if dt == 0.99*(dx^2)/(2.0*nu)
            fprintf('Diffusive CFL condition\n')
        elseif dt == 0.99*dx/max(abs(vx_new))
            fprintf('Convective CFL condition\n')
        end
        figure(1)
        semilogy(nxax,n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax) max(nxax)])
        hold on
        figure(2)
        plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(3)
        plot(vxax,vx_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(4)
        semilogy(nxax,n_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax) max(nxax)])
        hold on
        count = count + 1;
    end 
    
end

%%

figure(1)
xlabel('Position (m)','Fontsize',16)
ylabel('Density m$^{-3}$','Fontsize',16)
legend('show','Location','south')
hold off

figure(2)
xlabel('Position (m)','Fontsize',16)
ylabel('Mach number','Fontsize',16)
legend('show','Location','southeast')
hold off

figure(3)
xlabel('Position (m)','Fontsize',16)
ylabel('Velocity source ms$^{-1}$','Fontsize',16)
legend('show','Location','northwest')
hold off

figure(4)
xlabel('Position (m)','Fontsize',16)
ylabel('Density source m$^{-3}$','Fontsize',16)
legend('show','Location','south')
hold off




