%--------------------------------------------------------------------------------------------------------------%
% 1.5D solve for transport equations (continuity and momentum)                                                                           %
% NON-CONSERVATIVE FORMS                                                                                                                                                                                %
% staggered n and vz grids                                                                                     %
% momentum eqn upwind convection, central differenced diffusion (IMEX)                                         %
% continuity eqn first order upwind (explicit, flux selecting)                                                 %
% ghost points on density 
% nabla = (kapx, kapy, dz)
% rlbarnett c3149416 190723     
%--------------------------------------------------------------------------------------------------------------%

count = 2;
timerVal = tic;

input_transport_1_5D;

figure(1)
set(gcf,'Position',[563 925 560 420])
plot(nxax(2:npts-1),n_new(2:npts-1),'DisplayName',['time = 0s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density (m^{-3})','Fontsize',16)
legend('show','Location','west')
grid on
hold on

figure(2)
set(gcf,'Position',[7 925 560 420])
plot(vxax,vx_new/cs,'DisplayName',['time = 0s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Mach number','Fontsize',16)
legend('show','Location','southeast')
grid on
hold on
% 
% figure(3)
% set(gcf,'Position',[3 476 560 420])
% plot(vxax(2:npts-2),(vx_source(2:npts-2)),'DisplayName',['time = 0s'])
% xlabel('Position (m)','Fontsize',16)
% ylabel('Velocity source (ms^{-1})','Fontsize',16)
% legend('show','Location','northwest')
% grid on
% hold on
% 
% figure(4)
% set(gcf,'Position',[563 476 560 420])
% plot(nxax(2:npts-1),n_source(2:npts-1)*dt,'DisplayName',['time = 0s'])
% xlabel('Position (m)','Fontsize',16)
% ylabel('Density source (m^{-3})','Fontsize',16)
% legend('show','Location','northwest')
% grid on
% hold on

for ii=1:nmax

    vx = vx_new;
    
    rGhost = interp1([nxax(npts-2), nxax(npts-1)], [n_new(npts-2), n_new(npts-1)],...
        nxax(npts),'linear','extrap');   
    lGhost = interp1([nxax(2), nxax(3)], [n_new(2), n_new(3)],...
        nxax(1),'linear','extrap');
  
    for jj=2:npts-1
        if ((vx(1,jj-1)+vx(1,jj))/2)>0 
            nA(jj,jj) = - (1.0/ndx(1,jj-1))*vx(1,jj);
            nA(jj,jj-1) = (1.0/ndx(1,jj-1))*vx(1,jj-1);
        elseif ((vx(1,jj-1)+vx(1,jj))/2)<0
            nA(jj,jj) = (1.0/ndx(1,jj))*vx(1,jj-1);
            nA(jj,jj+1) = -(1.0/ndx(1,jj))*vx(1,jj);                
        end
        nA(jj,jj) = nA(jj,jj) + D_perp*(kap_x^2 + kap_y^2) -...
                (kap_x*vdrift_x + kap_y*vdrift_y);
    end

    An_exp = nI + dt*nA;

    An_exp(1,1) = 1.0;
    An_exp(end,end) = 1.0;

    n(1) = lGhost;
    n(1,end) = rGhost;

    n_new = An_exp*n' + dt*n_source';

    n_new = n_new';
    n = n_new;

    for jj=2:npts-2
        if vx(1,jj)>0
            vx_pos(jj,jj) = - (1.0/vdx(1,jj-1))*vx(1,jj) + (eta_perp/(mp*n(1,jj)) + D_perp)*...
                (kap_x^2 + kap_y^2) - (vdrift_x*kap_x + vdrift_y*kap*y) -... 
                (1.0/n(1,jj))*n_source(1,jj);
            vx_pos(jj,jj-1) = (1.0/vdx(1,jj-1))*vx(1,jj);
        elseif vx(1,jj)<0
            vx_neg(jj,jj) = (1.0/vdx(1,jj))*vx(1,jj) + (eta_perp/(mp*n(1,jj+1)) + D_perp)*...
                (kap_x^2 + kap_y^2) - (vdrift_x*kap_x + vdrift_y*kap*y) -... 
                (1.0/n(1,jj+1))*n_source(1,jj+1);
            vx_neg(jj,jj+1) = - (1.0/vdx(1,jj))*vx(1,jj);
        end
        vx_diff(jj,jj) = - (1.0/(vdx(1,jj-1)*vdx(1,jj)))*(2.0*eta_para);
        vx_diff(jj,jj-1) = (2.0/(vdx(1,jj-1)*(vdx(1,jj) + vdx(1,jj-1))))*eta_para;
        vx_diff(jj,jj+1) = (2.0/((vdx(1,jj-1) + vdx(1,jj))*vdx(1,jj)))*eta_para;
    end

    vxE = vx_pos + vx_neg;
    vxI = vx_diff;

    Avx_exp = vx_I + dt*vxE;
    Avx_imp = vx_I - dt*vxI;

    Avx_exp(1,1) = 1.0; Avx_exp(end,end) = 1.0;
    Avx_imp(1,1) = 1.0; Avx_imp(end,end) = 1.0;

    vx(1,1) = lvBC_val;
    vx(1,end) = rvBC_val;

    vx_source = source_stag(n,const.e,Te,Ti,const.mp,npts,ndx);
    pf_source = pond_source(const.me,m,om,const.e,Efield,vdx);
    pf_source = [0,pf_source,0];

    vx_source(1,1) = 0.0;
    vx_source(1,end) = 0.0;

    vx_newE = Avx_exp*vx';
    vx_new = Avx_imp\(vx_newE + dt*(vx_source'));

    vx_new = vx_new';

    nan_check = isnan(vx_new);
    
    if sum(nan_check) ~= 0
        fprintf('unstable, ii=%d\n',ii)
        return
    end

    if mod(ii,round(nmax/5))==0
        fprintf('***--------------------***\n')
        fprintf('ii=%d, count=%d\n', [ii count])
        fprintf('dt=%ds\n', dt)
        fprintf('total time=%ds\n', dt*ii)
        fprintf('simulation time %d\n', toc(timerVal))
        fprintf('Current rms tol calc %d\n', rms(vx - vx_new))
        figure(1)
        set(gcf,'Position',[563 925 560 420])
        plot(nxax,n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax+ndx(1)) max(nxax-ndx(end))])
        hold on
        figure(2)
        set(gcf,'Position',[7 925 560 420])
        plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(3)
        set(gcf,'Position',[3 476 560 420])
        plot(vxax(2:npts-2),(vx_source(2:npts-2)),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(4)
        set(gcf,'Position',[563 476 560 420])
        plot(nxax(2:npts-1),n_source(2:npts-1)*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlabel('Position (m)','Fontsize',16)
        ylabel('Density source ms^{-1}','Fontsize',16)
        legend('show','Location','northwest')
        hold on
        count = count + 1;
    end
    
end

fprintf('***--------------------***\n')
fprintf('ii=%d, count=%d\n', [ii count])
fprintf('dt=%ds\n', dt)
fprintf('total time=%ds\n', dt*ii)
fprintf('TOTAL simulation time %d\n', toc(timerVal))