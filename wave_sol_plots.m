% ----------------------plots----------------------- %

width = 1000;
height = 500;
x0 = 0;
y0 = 0;

figure(16)
set(gcf,'Position',[x0 y0 width height])
% set(gca,'YTickLabel','%.2f')

subplot(4,1,1)
plot(zax,source,'b','Linewidth',2)
ylabel('{\itJ_{y}} (Vm^{-1})')
set(gca, 'XTickLabel', [])
xlim([zmin,zmax])
set(gca,'Fontsize',20,'FontName','CMU Serif')
set(gcf,'Position',[x0 y0 width height],'Color','w')
ax = gca();
ax.YRuler.TickLabelFormat = '%1.1f';

subplot(4,1,2)
plot(zax, real(rf_ex), 'k','Linewidth',2)
ylabel('{\it E_x} (Vm^{-1})')

hold on

plot(zax, imag(rf_ex), '--r','Linewidth',2)
set(gca, 'XTickLabel', [])
% legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
xlim([zmin,zmax])
set(gca,'Fontsize',20,'FontName','CMU Serif')
set(gcf,'Position',[x0 y0+height width height],'Color','w')
ax = gca();
ax.YRuler.Exponent = 3;
ax.YRuler.TickLabelFormat = '%1.1f';

hold off

subplot(4,1,3)
plot(zax, real(rf_ey), 'k','Linewidth',2)
ylabel('{\it E_y} (Vm^{-1})')

hold on

plot(zax, imag(rf_ey), '--r','Linewidth',2)
set(gca, 'XTickLabel', [])
% legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
xlim([zmin,zmax]);
set(gca,'Fontsize',20,'FontName','CMU Serif')
set(gcf,'Position',[x0 y0+2*height width height],'Color','w')
ax = gca();
ax.YRuler.Exponent = 4;
ax.YRuler.TickLabelFormat = '%1.1f';

hold off

subplot(4,1,4)
plot(zax, real(rf_ez), 'k','Linewidth',2)
ylabel('{\it E_z} (Vm^{-1})')

hold on

plot(zax, imag(rf_ez), '--r','Linewidth',2)
xlabel('Position (m)')
% legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
xlim([zmin,zmax])
set(gca,'Fontsize',20,'FontName','CMU Serif')
set(gcf,'Position',[x0 y0+3*height width height],'Color','w')
ax = gca();
ax.YRuler.Exponent = 1;
ax.YRuler.TickLabelFormat = '%1.1f';

hold off

% export_fig('/Volumes/DATA/thesis/RFT/figs/vx_equilibrium.png',...
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

% export_fig('/Users/rhealbarnett/Documents/Documents/presentations/2020-rfscidac/abs_e_squ.png',...
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

