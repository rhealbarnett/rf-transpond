% ----------------------plots----------------------- %
x0 = 0;
y0 = 0;
width = 420*2;
height = 590*2;

figure(22)
set(gcf,'Units','points','Position',[x0 y0 width height],'Color','w')

subplot(7,1,1)
plot(zax,source,'k','linewidth',1.5)
xticks([])
set(gca,'Fontsize',18,'FontName','CMU Serif')
ylim([0, 2.2e6])
yticks([0.0, 1.0e6, 2.0e6])
yticklabels([{'0.0'} {'1.0'} {'2.0'}])
ylabel({'$J_y$'; '($\times 10^6$ Am$^{-2}$)'},'interpreter','latex')

subplot(7,1,2)
plot(zax, real(rf_ex), 'k','Linewidth',1.5)

hold on

plot(zax, imag(rf_ex), '--r','Linewidth',1.2)
plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-1.5e4, 1.5e4,npts),'color',...
    'm','linewidth',1)
plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-1.5e4, 1.5e4,npts),'color',...
    'm','linewidth',1)
set(gca, 'XTickLabel', [])
xlim([zmin zmax])
plot(zax,abs(rf_ex),'b','Linewidth',1.5)
set(gca,'Fontsize',18,'FontName','CMU Serif')
ax = gca();
ylim([-1.5e4, 1.5e4])
yticks([-1.5e4, 0.0, 1.5e4])
yticklabels([{'-1.5'} {'0.0'} {'1.5'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_x$'; '($\times 10^4$ Vm$^{-1}$)'},'interpreter','latex')


hold off

subplot(7,1,3)
plot(zax, real(rf_ey), 'k','Linewidth',1.5)

hold on

plot(zax, imag(rf_ey), '--r','Linewidth',1.2)
plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-6.0e4, 6.0e4,npts),'color',...
    'm','linewidth',1)
plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-6.0e4, 6.0e4,npts),'color',...
    'm','linewidth',1)
set(gca, 'XTickLabel', [])
xlim([zmin zmax])
plot(zax,abs(rf_ey),'b','Linewidth',1.5);
set(gca,'Fontsize',18,'FontName','CMU Serif')
ax = gca();
ylim([-6.0e4, 6.0e4])
yticks([-6.0e4, 0.0, 6.0e4])
yticklabels([{'-6.0'} {'0.0'} {'6.0'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_y$'; '($\times 10^4$ Vm$^{-1}$)'},'interpreter','latex')


hold off

subplot(7,1,4)
plot(zax, real(rf_ez), 'k','Linewidth',1.5)


hold on

plot(zax, imag(rf_ez), '--r','Linewidth',1.2)
plot(zax(floor(damp_len*npts))*ones(1,npts),linspace(-100,100,npts),'color',...
    'm','linewidth',1)
plot(-1*zax(floor(damp_len*npts))*ones(1,npts),linspace(-100,100,npts),'color',...
    'm','linewidth',1)
plot(zax,abs(rf_ez),'b','Linewidth',1.5)
xlim([zmin zmax])
set(gca,'Fontsize',18,'FontName','CMU Serif')
ax = gca();
xticks([])
ylim([-100,100])
yticks([-100,0.0,100])
yticklabels([{'-1.0'} {'0.0'} {'1.0'}])
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$E_z$'; '($\times 10^2$ Vm$^{-1}$)'},'interpreter','latex')

hold off

subplot(7,1,5)
plot(vxax,pf_source,'k','Linewidth',1.5)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',18,'FontName','CMU Serif')
xlim([zmin zmax])
ylim([-1.2e11,1.2e11])
yticks([-1.e11,0.0,1.e11])
yticklabels([{'-1.0'} {'0.0'} {'1.0'}])
ax = gca();
ax.YRuler.TickLabelFormat = '%1.1f';
ylabel({'$a_P$'; '($\times 10^{11}$ ms$^{-2}$)'},'interpreter','latex')

subplot(7,1,6)
plot(vxax,vx_new/cs,'k','Linewidth',1.5)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',18,'FontName','CMU Serif')
xlim([zmin zmax])
ax = gca();
ylabel({'Mach \#'},'interpreter','latex')
ax.YRuler.TickLabelFormat = '%1.1f';


subplot(7,1,7)
plot(nxax,n_new,'k','Linewidth',1.5)
set(gca,'Fontsize',18,'FontName','CMU Serif')
xlim([zmin zmax])
ylim([0.0, 1.5]*1e17)
yticks([0.0, 0.5e17, 1.0e17, 1.5e17])
yticklabels([{'0.0'} {'0.5'} {'1.0'} {'1.5'}])
ax = gca();
ylabel({'Density'; '($\times10^{17}$ m$^{-3}$)'},'interpreter','latex')
xlabel('Position (m)')

print('-f22',strcat(filepath,'outputs/7-panel-',num2str(ii)),'-dpng');
