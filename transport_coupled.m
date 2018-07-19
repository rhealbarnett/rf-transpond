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
transport_realistic;

%%

figure(1)
set(gcf,'Position',[563 925 560 420])
plot(nxax(2:npts-1),n_new,'DisplayName','time = 0s')
hold on

figure(2)
set(gcf,'Position',[7 925 560 420])
plot(vxax,vx_new/cs,'DisplayName','time = 0s')
hold on

figure(3)
set(gcf,'Position',[3 476 560 420])
plot(vxax(2:end-1),vx_source(2:end-1)*dt,'DisplayName','time = 0s')
hold on

figure(4)
set(gcf,'Position',[565 479 560 420])
plot(nxax(2:npts-1),n_source*dt,'DisplayName','time = 0s')
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
    
    alpha = dt/(dx);
    
    vx = vx_new;
    n = n_new;
%     nlog = log10(n);
    n_fit = interp1([nxax(2), nxax(3), nxax(npts-2), nxax(npts-1)],...
        [n(1), n(2), n(npts-3), n(npts-2)],...
        [nxax(1), nxax(npts)],'linear','extrap');
    n_extrap = [(n_fit(1)), n, (n_fit(2))];
    n_interp = interp1(nxax,n_extrap,vxax);
    gradn = (n_interp(3:end-1) - n_interp(2:end-2))/(dx);

    for jj=2:npts-2

%---------------------------------------------------------------------------%
% first order upwind                                                        %
% stable on staggered grid                                                  %
%---------------------------------------------------------------------------%
        if ((vx(1,jj-1)+vx(1,jj))/2)>0
            nA(jj,jj) = 1.0 - alpha*vx(1,jj);
            nA(jj,jj-1) = alpha*vx(1,jj-1);
            nA(npts-1,npts-1) = 1.0 - alpha*vx(1,end);
            nA(npts-1,npts-2) = alpha*vx(1,npts-2);
        elseif ((vx(1,jj-1)+vx(1,jj))/2)<0
            nA(jj,jj) = 1.0 + alpha*vx(1,jj-1);
            nA(jj,jj+1) = -alpha*vx(1,jj);
            nA(2,2) = 1.0 + alpha*vx(1,1);
            nA(2,3) = -alpha*vx(1,2);
        end 
        
%         n_source(jj-1,1) = n(1,jj-1)*n_neut(jj-1,1)*rate_coeff;
        
    end
    
%     n_source(1,1) = n_source(2,1);
%     n_source(end,1) = n_source(end-1,1);
    
    for jj=2:npts-2
        
        if vx(1,jj)>0
            vxA(jj,jj) = 1 - alpha*vx(1,jj);
            vxA(jj,jj-1) = alpha*vx(1,jj);
            vxA(end,end) = 1 - alpha*vx(1,end);
            vxA(end,end-1) = alpha*vx(1,end);
            pond_source(jj,1) = (1.0/m)*pond_const*((Efield(1,jj) - Efield(1,jj-1))/dx);
            vx_source(end,1) = -((Te + Ti)*e/(m*0.5*(n(1,end)+n(1,end-1))))*((n(1,end) - n(1,end-1))/dx) -...
            pond_source(jj,1);
        elseif vx(1,jj)<0
            vxA(jj,jj) = 1 + alpha*vx(1,jj);
            vxA(jj,jj+1) = -alpha*vx(1,jj);
            vxA(1,1) = 1 + alpha*vx(1,1);
            vxA(1,2) = -alpha*vx(1,1);
            pond_source(jj,1) = (1.0/m)*pond_const*((Efield(1,jj+1) - Efield(1,jj))/dx);
            vx_source(1,1) = -((Te + Ti)*e/(m*0.5*(n(1,1)+n(1,2))))*((n(1,2) - n(1,1))/dx) -...
            pond_source(jj,1);
        end
        
%         pond_source(jj,1) = (1.0/m)*pond_const*((Efield(1,jj+1) - Efield(1,jj-1)))/(2.0*dx);
%         vx_source(jj,1) = -((Te + Ti)*e)/(m*n_interp(1,jj))*(gradn(1,jj-1)) -...
%             pond_source(jj,1); 
%         pressure(1,jj) = (Te + Ti)*n_interp(1,jj-1)*e;
%         pressure_tot(1,jj) = pressure(1,jj) + (1/2)*n_interp(1,jj-1)*m*(vx(1,jj)^2);

        vx_source(jj,1) = -((Te + Ti)*e/(m*0.5*(n(1,jj)+n(1,jj-1))))*((n(1,jj) - n(1,jj-1))/dx) -...
    pond_source(jj,1);

    end
    
    nA = sparse(nA);
    vxA = sparse(vxA);
    
    n_new = dt*n_source + nA(2:npts-1,2:npts-1)*n';
    vx_new = dt*vx_source + vxA*vx';
    
    n_new = n_new';
    vx_new = vx_new';
     
    vx_new(1,1) = cs;
%     vx_new(1,end) = 0.0;
    
%     l_inf_vx(1,ii) = norm(vx - vx_new)/norm(vx);
%     l_two_vx(1,ii) = rms(vx - vx_new);
%     l_inf_n(1,ii) = norm(n - n_new)/norm(n);
%     l_two_n(1,ii) = rms(n - n_new);
%     
%     bound_check(1,ii) = gradn(end);
%     
    source_check(1,ii) = trapz(nxax(2:npts-1),n_source);
    
    flux = (vx_new.*n_interp);
    flux_check(ii,:) = flux;
% 
%     vx_mat(ii,:) = vx_new;
%     n_mat(ii,:) = n_new;
%     pressure_mat(ii,:) = pressure;
%     press_tot_mat(ii,:) = pressure_tot;
    
    nan_check = isnan(vx_new);
    
    if (0.008*(dx^2)/(2.0*nu))<(0.008*dx/max(abs(vx_new)))
        dt = 0.008*(dx^2)/(2.0*nu);
    elseif (0.008*(dx^2)/(2.0*nu))>(0.008*dx/max(abs(vx_new)))
        dt = 0.008*dx/max(abs(vx_new));
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
        if dt == 0.99*(dx^2)/(2.0*nu)
            fprintf('Diffusive CFL condition\n')
        elseif dt == 0.99*dx/max(abs(vx_new))
            fprintf('Convective CFL condition\n')
        end
        figure(1)
        set(gcf,'Position',[563 925 560 420])
        plot(nxax(2:end-1),n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax) max(nxax)])
        hold on
        figure(2)
        set(gcf,'Position',[7 925 560 420])
        plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(3)
        set(gcf,'Position',[3 476 560 420])
        plot(vxax(2:end-1),vx_source(2:end-1)*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(4)
        set(gcf,'Position',[565 479 560 420])
        plot(nxax(2:npts-1),n_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax) max(nxax)])
        hold on
        figure(5)
        set(gcf,'Position',[561 33 560 420])
        semilogy(tax(1:ii),source_check(1,1:ii))
        hold on
        semilogy(tax(1:ii),flux_check(1:ii,end),'r')
        semilogy(tax(1:ii),flux_check(1:ii,1),'*r')
        xlabel('Time (s)','Fontsize',16)
        ylabel('Particles m^{-2}s^{-1}','Fontsize',16)
        legend({'source','flux (right)','flux (left)'},'Fontsize',16)
        xlim([min(tax) max(tax)])
        hold off
        figure(6)
        set(gcf,'Position',[3 33 560 420])
        plot(vxax(2:npts-1),pond_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax) max(nxax)])
        hold on
        count = count + 1;
    end 
    
end

%%

figure(1)
set(gcf,'Position',[563 925 560 420])
plot(nxax(2:end-1),n_new,'DisplayName',['time = ' num2str(nmax*dt) ' s'])
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
plot(vxax(2:end-1),vx_source(2:end-1)*dt,'DisplayName',['time = ' num2str(nmax*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Velocity source ms^{-1}','Fontsize',16)
legend('show','Location','northwest')
hold off

figure(4)
set(gcf,'Position',[565 479 560 420])
plot(nxax(2:npts-1),n_source*dt,'DisplayName',['time = ' num2str(nmax*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density source m^{-3}','Fontsize',16)
legend('show','Location','south')
hold off

figure(5)
set(gcf,'Position',[561 33 560 420])
semilogy(tax,source_check(1,:))
hold on
semilogy(tax,flux_check(:,end),'r')
semilogy(tax,flux_check(:,1),'*r')
xlabel('Time (s)','Fontsize',16)
ylabel('Particles m^{-2}s^{-1}','Fontsize',16)
legend({'source','flux (right)','flux (left)'},'Fontsize',16)
xlim([min(tax) max(tax)])
hold off

figure(6)
set(gcf,'Position',[3 33 560 420])
plot(vxax(2:npts-1),pond_source*dt,'DisplayName',['time = ' num2str(nmax*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Ponderomotive source ms^{-1}','Fontsize',16)
legend('show','Location','south')
hold off

% figure(5)
% set(gcf,'Position',[561 33 560 420])
% semilogy(tax,l_inf_n);
% xlabel('Time (s)','Fontsize',16)
% ylabel('Relative difference in solution (for n)','Fontsize',16)
% 
% figure(6)
% set(gcf,'Position',[0 29 560 420])
% semilogy(tax,l_inf_vx);
% xlabel('Time (s)','Fontsize',16)
% ylabel('Relative difference in solution (for vx)','Fontsize',16)

% figure(7)
% levels = linspace(round(min(vx_mat(:)),-3),round(max(vx_mat(:)),-3),25);
% contourf(vxax,tax,vx_mat,levels,'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar
% 
% figure(8)
% levels = linspace(round(min(n_mat(:)),-3),round(max(n_mat(:)),-3),25);
% contourf(nxax,tax,n_mat,levels,'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar

% figure(9)
% levels = linspace((min(pressure_mat(:))),(max(pressure_mat(:))),25);
% contourf(vxax(1:npts-2),tax,pressure_mat,levels,'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar
% 
% figure(10)
% levels = linspace((min(press_tot_mat(:))),(max(press_tot_mat(:))),25);
% contourf(nxax(2:npts-1),tax,press_tot_mat,levels,'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar
