% ----------------------------------------------------------------------- %
% Script to read and plot matlab (.mat) files for the coupled 1D 
% wave / transport solver. 
%
% rlbarnett c3149416 191106
% ----------------------------------------------------------------------- %

data = 1;
equilibrium = 0;

if data
    
    filepath = '/Volumes/DATA/LAPD/matlab/results_jsource_ky20_n18_v2/';
%     load([filepath, 'coupled_transport_init.mat']);
    
% 
%     const = constants();
%     mach_init = equib.vx_new/equib.cs;
%     Te = equib.Te;
%     Ti = equib.Ti;
%     m = 4.00*const.amu;
%     e = const.e;
%     % cs = sqrt((Te+Ti)*e/m);
%     % mach_init = equib.vx_new./cs;
%     
%     Nmax = 1.0e17;
%     fact = Nmax/max(equib.n_new);
%     n_inter = fact*equib.n_new;
%     n_init = n_inter;
%     xax = linspace(equib.xmin,equib.xmax,equib.npts);

    lapd_equib;

    mach_init = vx_new./cs;
    
    width = 850;
    height = 1000;
    x0 = 0;
    y0 = 0;
    
    figure(1)
    set(gcf,'Position',[x0 y0 width height],'color','w')
    subplot(2,1,1)
    plot(nxax,n_new,'color',[0,0,0]+0.7,'Linewidth',5,...
        'DisplayName','time = 0T_{RF}')
    set(gca, 'XTickLabel', [])
    hold on
    box on

    subplot(2,1,2)
    plot(vxax,mach_init,'color',[0,0,0]+0.7,'Linewidth',5,...
        'DisplayName','time = 0T_{RF}')
    hold on
    box on
    
    figure(2)
    set(gcf,'Position',[x0 y0 width height],'color','w')
    subplot(3,1,1)
    plot(zax,real(rf_ex),'color',[0,0,0]+0.7,'Linewidth',5,...
        'DisplayName','time = 0T_{RF}')
    hold on
    plot(zax,imag(rf_ex),'color',[1,0,0,0.2],'Linewidth',5,...
        'DisplayName','time = 0T_{RF}')
    ax = gca();
    ax.YRuler.Exponent = 4;
    ax.YRuler.TickLabelFormat = '%1.1f';
    set(gca, 'XTickLabel', [])
    box on

    subplot(3,1,2)
    plot(zax,real(rf_ey),'color',[0,0,0]+0.7,'Linewidth',5,...
        'DisplayName','time = 0T_{RF}')
    hold on
    plot(zax,imag(rf_ey),'color',[1,0,0,0.2],'Linewidth',5,...
        'DisplayName','time = 0T_{RF}')
    ax = gca();
    ax.YRuler.Exponent = 4;
    ax.YRuler.TickLabelFormat = '%1.1f';
    set(gca, 'XTickLabel', [])
    box on

    subplot(3,1,3)
    plot(zax,real(rf_ez),'color',[0,0,0]+0.7,'Linewidth',5,...
        'DisplayName','time = 0T_{RF}')
    hold on
    plot(zax,imag(rf_ez),'color',[1,0,0,0.2],'Linewidth',5,...
        'DisplayName','time = 0T_{RF}')
    ax = gca();
    ax.YRuler.Exponent = 1;
    ax.YRuler.TickLabelFormat = '%1.1f';
    box on
    
    figure(3)
    subplot(2,1,1)
    plot(zax,abs(rf_ez).^2,'color',[0,0,0]+0.7,'Linewidth',5,'DisplayName',...
        ['time = 0T_{RF}'])
    hold on
    box on

%     figure(5)
    subplot(2,1,2)
    plot(vxax(2:npts-2),squeeze(pf_source(1,1,:)),'color',[0,0,0]+0.7,'Linewidth',5,'DisplayName',...
        ['time = 0T_{RF}'])
    hold on
    box on
    
    clear rf_ex rf_ey rf_ez rf_e
   

    for ii=nmax

        filename = strcat(filepath, 'coupled_transport_', num2str(ii),'.mat');

        load(filename)
        
        if exist('rf_e','var')
            rf_ez = rf_e(1,3:3:3*npts);
            rf_ey = rf_e(1,2:3:3*npts);
            rf_ex = rf_e(1,1:3:3*npts);
        else
        end

        figure(1)
        subplot(2,1,1)
        plot(nxax,n_new,'k','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        hold on
%         xlabel('Position (m)')
        ylabel('Density (m^{-3})')
        xlim([min(nxax) max(nxax)])
        ylim([min(n_new) max(n_new)+0.01*max(n_new)])
        set(gcf,'Position',[x0 y0 width height],'Color','w')
        set(gca,'Fontsize',30)
%         legend('show')
%         export_fig('/Volumes/DATA/LAPD/matlab/results_jsource_kyzero_v3/kyzero_v3_figs/sizetest.pdf','-r600')
%         print('/Volumes/DATA/LAPD/matlab/results_jsource_kyzero_n16_v2/n16_v2_figs/sizetest.pdf',...
%             '-dpdf','-r300')

%         figure(2)
        subplot(2,1,2)
        plot(vxax,vx_new/cs,'k','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        hold on
        xlabel('Position (m)')
        ylabel('Mach #')
        xlim([min(vxax) max(vxax)])
        ylim([-1 1])  
        set(gca,'Fontsize',30)
% 
%         export_fig('/Users/rhealbarnett/Documents/Documents/presentations/2020-rfscidac/results_nv_epara_ky20.png',...
%               '-r300')

%         legend('show')

        figure(2)
        subplot(3,1,1)
        plot(zax,real(rf_ex),'k','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        hold on
        plot(zax,imag(rf_ex),'--r','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
%         xlabel('Position (m)')
        ylabel('RF E_x (Vm^{-1})')
        xlim([min(zax) max(zax)])
        set(gca,'Fontsize',30)
%         legend('show')

        subplot(3,1,2)
        plot(zax,real(rf_ey),'k','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        hold on
        plot(zax,imag(rf_ey),'--r','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
%         xlabel('Position (m)')
        ylabel('RF E_y (Vm^{-1})')
        xlim([min(zax) max(zax)])
        set(gca,'Fontsize',30)

        subplot(3,1,3)
        plot(zax,real(rf_ez),'k','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        hold on
        plot(zax,imag(rf_ez),'--r','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        xlabel('Position (m)')
        ylabel('RF E_z (Vm^{-1})')
        xlim([min(zax) max(zax)])
        set(gca,'Fontsize',30)
%         legend('show')
%         export_fig('/Users/rhealbarnett/Documents/Documents/presentations/2020-rfscidac/results_efields_epara_ky20.png',...
%               '-r300')

        figure(3)
        subplot(2,1,1)
        plot(zax,abs(rf_ez).^2,'k','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        hold on
        xlabel('Position (m)')
        ylabel('|RF E_z (Vm^{-1})|^2')
        xlim([min(zax) max(zax)])
%         legend('show')

%         figure(5)
        subplot(2,1,2)
        plot(vxax,pf_source,'k','Linewidth',1.5,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        hold on
        xlabel('Position (m)')
        ylabel('PA (ms^{-2})')
        xlim([min(zax) max(zax)])
%         legend('show')



    end

elseif equilibrium
    
    filepath = '/Volumes/DATA/LAPD/matlab/equil_superrefined/';

    for ii=7066000:10000:7166000

        filename = strcat(filepath, 'equil_transport_', num2str(ii),'.mat');

        load(filename)
        
        figure(1)
        plot(nxax,n_new,'Linewidth',2,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        hold on
        xlabel('Position (m)')
        ylabel('Density (m^{-3})')
        xlim([min(nxax) max(nxax)])
        ylim([min(n_new) max(n_new)+0.01*max(n_new)])
        legend('show')

        figure(2)
        plot(vxax,vx_new/cs,'Linewidth',2,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        hold on
        xlabel('Position (m)')
        ylabel('Mach #')
        xlim([min(vxax) max(vxax)])
        ylim([-1 1])  
        legend('show')
        
        fprintf("Total number of particles in domain, %1.0f\n", trapz(n_new,nxax)) 
    
    end
    
end
