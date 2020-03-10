%% Calculate RF B amplitude from 1D cold plasma wave solver (wave_sol.m). 
%   Faradays Law : -i*om*B = Del X E
%   Bx = -(ky/om)*Ez - (i/om)dEy/dz
%   By = (i/om)*dEx/dz + (kx/om)*Ez
%   Bz = -(kx/om)*Ey + (ky/om)*Ex

function [rf_bx, rf_by, rf_bz] = B_RF(ax,kx,ky,om,rf_ex,rf_ey,rf_ez,plots,source)

    npts = length(ax);
    axmax = ax(1,end);
    axmin = ax(1,1);
    h = (axmax - axmin)/(npts-1);

    rfex_diff = zeros(1,npts);
    rfey_diff = zeros(1,npts);

    rfex_diff(2:npts-1) = (rf_ex(3:npts) - rf_ex(1:npts-2)) / (2.0*h);
    rfex_diff(1) = (-3.0*rf_ex(1) + 4.0*rf_ex(2) - rf_ex(3)) / (2.0*h);
    rfex_diff(npts) = (rf_ex(npts-2) - 4.0*rf_ex(npts-1) + 3.0*rf_ex(npts)) / (2.0*h);

    rfey_diff(2:npts-1) = (rf_ey(3:npts) - rf_ey(1:npts-2)) / (2.0*h);
    rfey_diff(1) = (-3.0*rf_ey(1) + 4.0*rf_ey(2) - rf_ey(3)) / (2.0*h);
    rfey_diff(npts) = (rf_ey(npts-2) - 4.0*rf_ey(npts-1) + 3.0*rf_ey(npts)) / (2.0*h);

    rf_bx = (ky/om).*rf_ez + (1i/om)*rfey_diff;
    rf_by = -(1i/om)*rfex_diff - (kx./om).*rf_ez;
    rf_bz = (kx./om).*rf_ey - (ky./om).*rf_ex;
    
    if plots
        
        width = 1000;
        height = 1000;
        x0 = 0;
        y0 = 0;
        
        figure(13)
        set(gcf,'Position',[x0 y0 width height],'color','w')
        subplot(3,2,2)
        plot(ax, real(rf_bx), 'k','Linewidth',2)
        ylabel('{\it B_x} (T)')

        hold on

        plot(ax, imag(rf_bx), '--r','Linewidth',2)
        set(gca, 'XTickLabel', [])
%         legend('Re[Bx]', 'Im[Bx]', 'Location', 'northwest')
        xlim([axmin,axmax])
        set(gca,'Fontsize',20)
        axes = gca();
%         axes.YRuler.Exponent = -8;
%         axes.YRuler.TickLabelFormat = '%1.1f';
        

        hold off

        subplot(3,2,4)
        plot(ax, real(rf_by), 'k','Linewidth',2)
        ylabel('{\it B_y} (T)')

        hold on

        plot(ax, imag(rf_by), '--r','Linewidth',2)
        set(gca, 'XTickLabel', [])
%         legend('Re[By]', 'Im[By]', 'Location', 'northwest')
        xlim([axmin,axmax]);
        set(gca,'Fontsize',20)
        axes = gca();
        
%         axes.YRuler.Exponent = -10;
%         axes.YRuler.TickLabelFormat = '%1.1f';

        hold off

        subplot(3,2,6)
        plot(ax, real(rf_bz), 'k','Linewidth',2)
        ylabel('{\it B_z} (T)')

        hold on

        plot(ax, imag(rf_bz), '--r','Linewidth',2)
        xlabel('Position (m)')
%         legend('Re[Bz]', 'Im[Bz]', 'Location', 'northwest')
        xlim([axmin,axmax])
        set(gca,'Fontsize',20)
        axes = gca();
%         axes.YRuler.Exponent = -8;
%         axes.YRuler.TickLabelFormat = '%1.1f';

        hold off
        
        subplot(3,2,1)
        plot(ax, real(rf_ex), 'k','Linewidth',2)
        ylabel('{\it E_x} (Vm^{-1})')

        hold on

        plot(ax, imag(rf_ex), '--r','Linewidth',2)
        set(gca, 'XTickLabel', [])
%         legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
        xlim([axmin,axmax])
        set(gca,'Fontsize',20)
        axes = gca();
%         axes.YRuler.Exponent = -3;
%         axes.YRuler.TickLabelFormat = '%1.1f';

        hold off

        subplot(3,2,3)
        plot(ax, real(rf_ey), 'k','Linewidth',2)
        ylabel('{\it E_y} (Vm^{-1})')

        hold on

        plot(ax, imag(rf_ey), '--r','Linewidth',2)
        set(gca, 'XTickLabel', [])
%         legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
        xlim([axmin,axmax]);
        set(gca,'Fontsize',20)
        axes = gca();
%         axes.YRuler.Exponent = -2;
%         axes.YRuler.TickLabelFormat = '%1.1f';

        hold off

        subplot(3,2,5)
        plot(ax, real(rf_ez), 'k','Linewidth',2)
        ylabel('{\it E_z} (Vm^{-1})')

        hold on

        plot(ax, imag(rf_ez), '--r','Linewidth',2)
        xlabel('Position (m)')
%         legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
        xlim([axmin,axmax])
        set(gca,'Fontsize',20)
        axes = gca();
%         axes.YRuler.Exponent = -5;
%         axes.YRuler.TickLabelFormat = '%1.1f';

        hold off
        
        width = 600;
        height = 700;
        
        figure(14)
        set(gcf,'Position',[x0 y0 width height],'color','w')
        
        subplot(4,1,1)
        plot(ax,source,'b','Linewidth',2)
        ylabel('{\iti}\omega\mu_0{\it J_{y}} (Vm^{-1})')
        set(gca, 'XTickLabel', [])
        set(gca,'Fontsize',22,'FontName','CMU Serif')
%         ax = gca();
%         ax.YRuler.TickLabelFormat = '%1.1f';
        
        subplot(4,1,2)
        
        plot(ax, real(rf_bx), 'k','Linewidth',2)
        ylabel('{\it B_x} (T)')

        hold on

        plot(ax, imag(rf_bx), '--r','Linewidth',2)
        set(gca, 'XTickLabel', [])
%         legend('Re[Bx]', 'Im[Bx]', 'Location', 'northwest')
        xlim([axmin,axmax])
        set(gca,'Fontsize',22)
        axes = gca();
%         plot(-0.06*ones(1,npts),linspace(-5e-8,5e-8,npts),...
%             '--b','Linewidth',1.5)
%         plot(0.06*ones(1,npts),linspace(-5e-8,5e-8,npts),...
%             '--b','Linewidth',1.5)
%         axes.YRuler.Exponent = -8;
%         axes.YRuler.TickLabelFormat = '%1.1f';

        hold off

        subplot(4,1,3)
        plot(ax, real(rf_by), 'k','Linewidth',2)
        ylabel('{\it B_y} (T)')

        hold on

        plot(ax, imag(rf_by), '--r','Linewidth',2)
        set(gca, 'XTickLabel', [])
%         legend('Re[By]', 'Im[By]', 'Location', 'northwest')
        xlim([axmin,axmax]);
        set(gca,'Fontsize',22)
        axes = gca();
%         plot(-0.06*ones(1,npts),linspace(-6e-10,6e-10,npts),...
%             '--b','Linewidth',1.5)
%         plot(0.06*ones(1,npts),linspace(-6e-10,6e-10,npts),...
%             '--b','Linewidth',1.5)
%         axes.YRuler.Exponent = -10;
%         axes.YRuler.TickLabelFormat = '%1.1f';

        hold off

        subplot(4,1,4)
        plot(ax, real(rf_bz), 'k','Linewidth',2)
        ylabel('{\it B_z} (T)')

        hold on

        plot(ax, imag(rf_bz), '--r','Linewidth',2)
        xlabel('Position (m)')
%         legend('Re[Bz]', 'Im[Bz]', 'Location', 'northwest')
        xlim([axmin,axmax])
        set(gca,'Fontsize',22)
        axes = gca();
%         plot(-0.06*ones(1,npts),linspace(-5e-8,5e-8,npts),...
%             '--b','Linewidth',1.5)
%         plot(0.06*ones(1,npts),linspace(-5e-8,5e-8,npts),...
%             '--b','Linewidth',1.5)
%         axes.YRuler.Exponent = -8;
%         axes.YRuler.TickLabelFormat = '%1.1f';

        hold off
        
%         export_fig('/Users/rhealbarnett/Documents/Documents/presentations/2020-rfscidac/b_rf_wsource_zoomed.png',...
%             '-r300')
        
    end

end
