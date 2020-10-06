% ----------------------plots----------------------- %

x0 = 0;
y0 = 0;
width = 420*2;
height = 590*2;

col_purp = [0.28 0.06 0.47]; 

figure(17)
set(gcf,'Position',[x0 y0 width height])
% set(gca,'YTickLabel','%.2f')

subplot(4,1,1)
plot(zax,source,'b','Linewidth',1.5)
ylabel('{\itJ_{y}} (Vm^{-1})')
set(gca, 'XTickLabel', [])
xlim([zmin,zmax])
set(gca,'Fontsize',20,'FontName','CMU Serif')
set(gcf,'Position',[x0 y0 width height],'Color','w')
ax = gca();
% ax.YRuler.TickLabelFormat = '%1.1f';

subplot(4,1,2)
plot(zax, real(rf_ex), 'k','Linewidth',1.5)
ylabel('{\it E_x} (Vm^{-1})')

hold on

plot(zax, imag(rf_ex), '--r','Linewidth',1.5)
plot(zax, abs(rf_ex),'b','Linewidth',1.5)
set(gca, 'XTickLabel', [])
% legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
xlim([zmin,zmax])
set(gca,'Fontsize',20,'FontName','CMU Serif')
set(gcf,'Position',[x0 y0+height width height],'Color','w')
ax = gca();
% ax.YRuler.Exponent = 3;
% ax.YRuler.TickLabelFormat = '%1.1f';

hold off

subplot(4,1,3)
plot(zax, real(rf_ey), 'k','Linewidth',1.5)
ylabel('{\it E_y} (Vm^{-1})')

hold on

plot(zax, imag(rf_ey), '--r','Linewidth',1.5)
plot(zax, abs(rf_ey), 'b','Linewidth',1.5)
set(gca, 'XTickLabel', [])
% legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
xlim([zmin,zmax]);
% ylim([-0.2 0.2])
set(gca,'Fontsize',20,'FontName','CMU Serif')
set(gcf,'Position',[x0 y0+2*height width height],'Color','w')
ax = gca();
% ax.YRuler.Exponent = 4;
% ax.YRuler.TickLabelFormat = '%1.1f';

hold off

subplot(4,1,4)
plot(zax, real(rf_ez), 'k','Linewidth',1.5)
ylabel('{\it E_z} (Vm^{-1})')

hold on

plot(zax, imag(rf_ez), '--r','Linewidth',1.5)
plot(zax, abs(rf_ez), 'b','Linewidth',1.5)
xlabel('Position (m)')
% ylim([-0.65,0.65])
% legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
xlim([zmin,zmax])
set(gca,'Fontsize',20,'FontName','CMU Serif')
set(gcf,'Position',[x0 y0+3*height width height],'Color','w')
ax = gca();
% ax.YRuler.Exponent = 1;
% ax.YRuler.TickLabelFormat = '%1.1f';

hold off

% export_fig('/Volumes/DATA/matlab/wave_verification/smithe_email/efields_jxjy_polcalc_decay.png',...
%     '-r300')


%%

figure(17)

x0 = 0;
y0 = 0;
width = 700;
height = 500;

subplot(3,1,1)
set(gcf,'Position',[x0 y0 width height],'color','w')
plot(zax, abs(rf_ex).^2, 'k','Linewidth',2)
ylabel('|E_x|^2 (V^2m^{-2})')
set(gca, 'XTickLabel', [])
legend('Ex', 'Location', 'northwest')
xlim([zmin,zmax])
set(gca,'Fontsize',20)


subplot(3,1,2)
plot(zax, abs(rf_ey).^2, 'k','Linewidth',2)
ylabel('|E_y|^2 (V^2m^{-2})')
set(gca, 'XTickLabel', [])
legend('Ey', 'Location', 'northwest')
xlim([zmin,zmax]);
set(gca,'Fontsize',20)


subplot(3,1,3)
plot(zax, abs(rf_ez).^2, 'k','Linewidth',2)
ylabel('|E_z|^2 (V^2m^{-2})')
legend('Ez', 'Location', 'northwest')
xlim([zmin,zmax]);
set(gca,'Fontsize',20)

% export_fig('/Volumes/DATA/matlab/scan_results/figs/gaussian_eperp_100kW_15i_2e17_damplines.png',...
%     '-r300')

%%

pa_epara = squeeze(pf(1,1,:));
pa_ipara = squeeze(pf(1,2,:));
pa_eperp = squeeze(pf(2,1,:));
pa_iperp = squeeze(pf(2,2,:));

figure(18)
set(gcf,'Position',[x0 y0 width height],'color','w')
subplot(2,2,1)
plot(zax, pa_epara, 'k','Linewidth',2)
ylabel('PA_{e,||} (ms^{-2})')
set(gca, 'XTickLabel', [])
% legend('Ex', 'Location', 'northwest')
xlim([zmin,zmax])
set(gca,'Fontsize',20)


subplot(2,2,2)
plot(zax, pa_ipara, 'k','Linewidth',2)
ylabel('PA_{He,||} (ms^{-2})')
set(gca, 'XTickLabel', [])
% legend('Ey', 'Location', 'northwest')
xlim([zmin,zmax]);
set(gca,'Fontsize',20)


subplot(2,2,3)
plot(zax, pa_eperp, 'k','Linewidth',2)
ylabel('PA_{e,\perp} (ms^{-2})')
xlabel('Position (m)')
xlim([zmin,zmax]);
set(gca,'Fontsize',20)

subplot(2,2,4)
plot(zax, pa_iperp, 'k','Linewidth',2)
ylabel('PA_{He,\perp} (ms^{-2})')
xlabel('Position (m)')
xlim([zmin,zmax]);
set(gca,'Fontsize',20)

% export_fig('/Users/rhealbarnett/Documents/Documents/presentations/2020-rfscidac/pond_accel.png',...
%     '-r300')

%%

x0 = 0;
y0 = 0;
width = 420*2;
height = 590*2;

col_purp = [0.28 0.06 0.47];            
col_pink = [0.78 0.24 0.3];
col_oran = [0.94 0.44 0.13];

figure(21)
set(gcf,'Units','points','Position',[x0 y0 width height],'Color','w')
% set(gca,'YTickLabel','%.2f')

subplot(7,1,1)
plot(zax,source,'k','Linewidth',1.5)
set(gca, 'XTickLabel', [])
xlim([zmin zmax])
set(gca,'Fontsize',18,'FontName','CMU Serif')
% set(gcf,'Position',[x0 y0 width height],'Color','w')
ax = gca();
ax.YRuler.TickLabelFormat = '%1.1f';
% ax.YRuler.Exponent = 5;
ylim([0.0 2.7e6])
yticks([0.0 2.7e6])
yticklabels([{'0.0'} {'2.0'}])
ylabel({'$J_y$'; '($\times 10^6$ Am$^{-2}$)'},'interpreter','latex')
text(0.005,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black')

subplot(7,1,2)
plot(zax, real(rf_ex), 'k','Linewidth',1.5)

hold on

plot(zax, imag(rf_ex), '--r','Linewidth',1.2)
plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-5.0e2,5.0e2,npts),'color',...
    col_purp,'linewidth',1)
plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-5.0e2,5.0e2,npts),'color',...
    col_purp,'linewidth',1)
set(gca, 'XTickLabel', [])
% legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
xlim([zmin zmax])
plot(zax,abs(rf_ex),'b','Linewidth',1.5)
set(gca,'Fontsize',18,'FontName','CMU Serif')
% set(gcf,'Position',[x0 y0+height width height],'Color','w')
ax = gca();
% % ax.YRuler.Exponent = 1;
ylim([-5.2e2, 5.2e2])
yticks([-5.0e2, 0.0, 5.0e2])
yticklabels([{'-5.0'} {'0.0'} {'5.0'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_x$'; '($\times 10^2$ Vm$^{-1}$)'},'interpreter','latex')
text(0.005,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black')

hold off

subplot(7,1,3)
plot(zax, real(rf_ey), 'k','Linewidth',1.5)

hold on

plot(zax, imag(rf_ey), '--r','Linewidth',1.2)
plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-2.5e3,2.5e3,npts),'color',...
    col_purp,'linewidth',1)
plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-2.5e3,2.5e3,npts),'color',...
    col_purp,'linewidth',1)
set(gca, 'XTickLabel', [])
% legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
xlim([zmin zmax])
plot(zax,abs(rf_ey),'b','Linewidth',1.5);
set(gca,'Fontsize',18,'FontName','CMU Serif')
% set(gcf,'Position',[x0 y0+2*height width height],'Color','w')
ax = gca();
% % ax.YRuler.Exponent = 3;
ylim([-2.5e3, 2.5e3])
yticks([-2.5e3, 0.0, 2.5e3])
yticklabels([{'-2.0'} {'0.0'} {'2.0'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_y$'; '($\times 10^3$ Vm$^{-1}$)'},'interpreter','latex')
text(0.005,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black')

hold off

subplot(7,1,4)
plot(zax, real(rf_ez), 'k','Linewidth',1.5)


hold on

plot(zax, imag(rf_ez), '--r','Linewidth',1.2)
plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-5,5,npts),'color',...
    col_purp,'linewidth',1)
plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-5,5,npts),'color',...
    col_purp,'linewidth',1)
% xlabel('Position (m)')
plot(zax,abs(rf_ez),'b','Linewidth',1.5)
% legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
xlim([zmin zmax])
set(gca,'Fontsize',18,'FontName','CMU Serif')
% set(gcf,'Position',[x0 y0+3*height width height],'Color','w')
ax = gca();
set(gca, 'XTickLabel', [])
% % ax.YRuler.Exponent = 0;
ylim([-5,5])
yticks([-5,0.0,5])
yticklabels([{'-5.0'} {'0.0'} {'5.0'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_z$'; '(Vm$^{-1}$)'},'interpreter','latex')
text(0.005,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black')

hold off
% 
% 
% pa_epara = squeeze(pond(1,1,:));
% pa_ipara = squeeze(pond(1,2,:));
% pa_eperp = squeeze(pond(2,1,:));
% pa_iperp = squeeze(pond(2,2,:));

subplot(7,1,5)
plot(vxax,pf_source,'k','Linewidth',1.5)
% plot(zax,(pa_epara+pa_ipara+pa_eperp+pa_iperp),'k','Linewidth',2)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',18,'FontName','CMU Serif')
xlim([zmin zmax])
ylim([-1.e6,1.e6])
yticks([-1.e6,0.0,1.e6])
yticklabels([{'-1.0'} {'0.0'} {'1.0'}])
ax = gca();
% ax.YRuler.Exponent = 8;
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$a_P$'; '($\times 10^6$ ms$^{-2}$)'},'interpreter','latex')
text(0.005,0.98,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black')

subplot(7,1,6)
plot(vxax,vx_new/cs,'k','Linewidth',1.5)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',18,'FontName','CMU Serif')
xlim([zmin zmax])
ax = gca();
% ax.YRuler.Exponent = 0;
ylabel({'Mach \#'},'interpreter','latex')
ax.YRuler.TickLabelFormat = '%1.1f';
% text(0.005,0.98,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
%                 'color','black')

subplot(7,1,7)
plot(nxax,(n_init - n_new)./n_init,'k','Linewidth',1.5)
set(gca,'Fontsize',18,'FontName','CMU Serif')
xlim([zmin zmax])
ylim([-4., 6.]*1e-5)
yticks([-4.e-5, 0.0, 4.e-5])
yticklabels([{'-4.0'} {'0.0'} {'4.0'}])
ax = gca();
% ax.YRuler.Exponent = -2;
% ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$R_n \times10^{-5}$'},'interpreter','latex')
xlabel('Position (m)')
text(0.005,0.98,'(g)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black')

% export_fig('/Volumes/DATA/thesis/figs/delta_onespike_epara_100kW_n1e17.png',...
%     '-r300')

%%

x0 = 0;
y0 = 0;
width = 1500;
height = 700;

figure(22)
set(gcf,'Units','points','Position',[x0 y0 width height],'Color','w')
% set(gca,'YTickLabel','%.2f')

subplot(3,2,1)
plot(zax, real(rf_ex), 'k','Linewidth',1.5)

hold on

plot(zax, imag(rf_ex), '--r','Linewidth',1.2)
% plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-5.0e2,5.0e2,npts),'color',...
%     col_purp,'linewidth',1)
% plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-5.0e2,5.0e2,npts),'color',...
%     col_purp,'linewidth',1)
set(gca, 'XTickLabel', [])
% legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
xlim([zmin zmax])
plot(zax,abs(rf_ex),'b','Linewidth',1.5)
set(gca,'Fontsize',25,'FontName','CMU Serif')
% set(gcf,'Position',[x0 y0+height width height],'Color','w')
ax = gca();
% % ax.YRuler.Exponent = 1;
ylim([-5.2e2, 5.2e2])
yticks([-5.0e2, 0.0, 5.0e2])
yticklabels([{'-5.0'} {'0.0'} {'5.0'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_x$'; '($\times 10^2$ Vm$^{-1}$)'},'interpreter','latex')


hold off

subplot(3,2,3)
plot(zax, real(rf_ey), 'k','Linewidth',1.5)

hold on

plot(zax, imag(rf_ey), '--r','Linewidth',1.2)
% plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-2.5e3,2.5e3,npts),'color',...
%     col_purp,'linewidth',1)
% plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-2.5e3,2.5e3,npts),'color',...
%     col_purp,'linewidth',1)
set(gca, 'XTickLabel', [])
% legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
xlim([zmin zmax])
plot(zax,abs(rf_ey),'b','Linewidth',1.5);
set(gca,'Fontsize',25,'FontName','CMU Serif')
% set(gcf,'Position',[x0 y0+2*height width height],'Color','w')
ax = gca();
% % ax.YRuler.Exponent = 3;
ylim([-2.5e3, 2.5e3])
yticks([-2.5e3, 0.0, 2.5e3])
yticklabels([{'-2.0'} {'0.0'} {'2.0'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_y$'; '($\times 10^3$ Vm$^{-1}$)'},'interpreter','latex')


hold off

subplot(3,2,5)
plot(zax, real(rf_ez), 'k','Linewidth',1.5)


hold on

plot(zax, imag(rf_ez), '--r','Linewidth',1.2)
% plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-5,5,npts),'color',...
%     col_purp,'linewidth',1)
% plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-5,5,npts),'color',...
%     col_purp,'linewidth',1)
% xlabel('Position (m)')
plot(zax,abs(rf_ez),'b','Linewidth',1.5)
% legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
xlim([zmin zmax])
set(gca,'Fontsize',25,'FontName','CMU Serif')
% set(gcf,'Position',[x0 y0+3*height width height],'Color','w')
ax = gca();
% set(gca, 'XTickLabel', [])
xlabel('Position (m)')
% % ax.YRuler.Exponent = 0;
ylim([-5,5])
yticks([-5,0.0,5])
yticklabels([{'-5.0'} {'0.0'} {'5.0'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_z$'; '(Vm$^{-1}$)'},'interpreter','latex')

hold off

subplot(3,2,2)
plot(vxax,pf_source,'k','Linewidth',1.5)
% plot(zax,(pa_epara+pa_ipara+pa_eperp+pa_iperp),'k','Linewidth',2)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',25,'FontName','CMU Serif')
xlim([zmin zmax])
ylim([-1.2e8,1.2e8])
yticks([-1.e8,0.0,1.e8])
yticklabels([{'-1.0'} {'0.0'} {'1.0'}])
ax = gca();
% ax.YRuler.Exponent = 8;
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$a_P$'; '($\times 10^8$ ms$^{-2}$)'},'interpreter','latex')

subplot(3,2,4)
plot(vxax,vx_new/cs,'k','Linewidth',1.5)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',25,'FontName','CMU Serif')
xlim([zmin zmax])
ax = gca();
% ax.YRuler.Exponent = 0;
ylabel({'Mach \#'},'interpreter','latex')
ax.YRuler.TickLabelFormat = '%1.1f';


subplot(3,2,6)
plot(nxax,(n_init - n_new)./n_init,'k','Linewidth',1.5)
set(gca,'Fontsize',25,'FontName','CMU Serif')
xlim([zmin zmax])
ylim([-1.5, 1.5]*1e-2)
yticks([-1.0e-2, 0.0, 1.0e-2])
yticklabels([{'-1.0'} {'0.0'} {'1.0'}])
ax = gca();
% ax.YRuler.Exponent = -2;
% ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$R_n \times10^{-2}$'},'interpreter','latex')
xlabel('Position (m)')
