% ----------------------------------------------------------------------- %
% Script to read and plot matlab (.mat) files for the coupled 1D 
% wave / transport solver. 
%
% rlbarnett c3149416 191106
% ----------------------------------------------------------------------- %

data = 1;
equilibrium = 0;

width = 600;
height = 400;
x0 = 0;
y0 = 0;

if data
    
    equib = load('/Volumes/DATA/LAPD/matlab/inputs/equil_transport_input.mat','vx_new','n_new','xmin',...
        'xmax','npts','Te','Ti','nxax','vxax','cs');

    const = constants();
    mach_init = equib.vx_new/equib.cs;
    Te = equib.Te;
    Ti = equib.Ti;
    m = 4.00*const.amu;
    e = const.e;
    % cs = sqrt((Te+Ti)*e/m);
    % mach_init = equib.vx_new./cs;
    n_init = equib.n_new;
    xax = linspace(equib.xmin,equib.xmax,equib.npts);
    
    figure(1)
    plot(equib.nxax,n_init,'color',[0,0,0]+0.8,'Linewidth',4,...
        'DisplayName','time = 0s')
    hold on
    box on

    figure(2)
    plot(equib.vxax,mach_init,'color',[0,0,0]+0.8,'Linewidth',4,...
        'DisplayName','time = 0s')
    hold on
    box on
    
    clear equib
    
    filepath = '/Volumes/DATA/LAPD/matlab/results_jsource_kyzero_v3/';

    for ii=106

        filename = strcat(filepath, 'coupled_transport_', num2str(ii),'.mat');

        load(filename)
        
        if exist('rf_e','var')
            rf_ez = rf_e(1,3:3:3*npts);
        else
        end

        figure(1)
        plot(nxax,n_new,'k','Linewidth',2,'DisplayName',...
            ['time = ' num2str(round(double(ii)*dt/period)) 'T_{RF}'])
        hold on
        xlabel('Position (m)','Fontsize',20,'FontName','CMU Serif')
        ylabel('Density (m^{-3})')
        xlim([min(nxax) max(nxax)])
        ylim([min(n_new) max(n_new)+0.01*max(n_new)])
        set(gcf,'Position',[x0 y0 width height],'Color','w')
        legend('show')
%         export_fig('/Volumes/DATA/LAPD/matlab/results_jsource_kyzero_v3/kyzero_v3_figs/sizetest.pdf','-r600')
%         print('/Volumes/DATA/LAPD/matlab/results_jsource_kyzero_n16_v2/n16_v2_figs/sizetest.pdf',...
%             '-dpdf','-r300')

        figure(2)
        plot(vxax,vx_new/cs,'k','Linewidth',2,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        hold on
        xlabel('Position (m)')
        ylabel('Mach #')
        xlim([min(vxax) max(vxax)])
        ylim([-1 1])  
        legend('show')

        figure(3)
        plot(xax,real(rf_ez),'k','Linewidth',2,'DisplayName',['Re[Ez] time = ' num2str(double(ii)*dt) ' s'])
        hold on
        plot(xax,imag(rf_ez),'--r','Linewidth',2,'DisplayName',['Im[Ez] time = ' num2str(double(ii)*dt) ' s'])
        xlabel('Position (m)')
        ylabel('RF E_z (Vm^{-1})')
        xlim([min(xax) max(xax)])
        legend('show')

        figure(4)
        plot(xax,abs(rf_ez).^2,'k','Linewidth',2,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        hold on
        xlabel('Position (m)')
        ylabel('|RF E_z (Vm^{-1})|^2')
        xlim([min(xax) max(xax)])
        legend('show')

        figure(5)
        plot(vxax,pf_source,'k','Linewidth',2,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        hold on
        xlabel('Position (m)')
        ylabel('PF source (ms^{-2})')
        xlim([min(xax) max(xax)])
        legend('show')

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
