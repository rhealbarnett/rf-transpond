%-------------------------------------------------------------------------%
% Projection plots. 
% rlbarnett c3149416 200109
%-------------------------------------------------------------------------%

width = 1600;
height = 500;
x0 = 0;
y0 = 0;

figure(1)
ax1 = subplot(1,3,1);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-2.5e-3,2.5e-3,50);
contourf(xax,zax,(Bx_x),levels,'Linecolor','none')
c = colorbar;
c.Ticks = ([-2e-3, -1.e-3, 0.0, 1.0e-3, 2.0e-3]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
% map = redblue(256);
% colormap(ax1,map(1:256/2,:))
colormap(ax1,redblue)
ylim([-0.1 0.1])
% caxis([-250e1 250e1])
% xlabel('z axis Position (m)')
ylabel('$z$ Position (m)','interpreter','latex','Fontsize',30)
xlabel('$x$ Position (m)','interpreter','latex','Fontsize',30)
ylabel(c,'$B_x(x,z)$ ($\times 10^{-3}$ Vm$^{-1}$)','interpreter','latex','Fontsize',30)
xlim([0. 0.215])
set(gca, 'XDir','reverse')
set(gca,'Linewidth',2)
% ax = gca();
% c.Ruler.Exponent = 2;
% c.Ruler.TickLabelFormat = '%1.f';


ax2 = subplot(1,3,2);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-2.5e-3,2.5e-3,50);
contourf(xax,zax,(By_x),levels,'Linecolor','none')
colormap(ax2,redblue)
caxis([-2.5e-3 2.4e-3])
% xlabel('z axis Position (m)')
xlabel('$x$ Position (m)','interpreter','latex','Fontsize',30)
ylim([-0.1 0.1])
c = colorbar;
c.Ticks = ([-2e-5, -1.e-5, 0.0, 1.0e-5, 2.0e-5]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
ylabel(c,'$B_y(x,z)$ ($\times 10^{-5}$ Vm$^{-1}$)','interpreter','latex','Fontsize',30)
% ax = gca();
xlim([0. 0.215])
set(gca, 'XDir','reverse')
% caxis([-1.1e4 1.1e4])
% c.Ruler.Exponent = 2;
% c.Ruler.TickLabelFormat = '%1.f';
set(gca,'Linewidth',2)
set(gca,'ytick',[])

ax3 = subplot(1,3,3);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-2.5e-3,2.5e-3,50);
contourf(xax,zax,(Bz_x),levels,'Linecolor','none')
colormap(ax3,redblue)
xlabel('$x$ Position (m)','interpreter','latex','Fontsize',30)
% ylabel('x Position (m)')
c = colorbar;
c.Ticks = ([-2e-3, -1.e-3, 0.0, 1.0e-3, 2.0e-3]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
ylim([-0.1 0.1])
ylabel(c,'$B_z(x,z)$ ($\times 10^{-3}$ Vm$^{-1}$)','interpreter','latex','Fontsize',30)
xlim([0. 0.215])
set(gca,'ytick',[])
set(gca, 'XDir','reverse')
set(gca,'Linewidth',2)
% c.Ruler.Exponent = 3;
% caxis([-2.0e2 2.0e2])
% c.Ruler.TickLabelFormat = '%1.f';

% export_fig('/Volumes/DATA/LAPD/matlab/wave_projection/e_xz_negxax.png',...
%     '-r300')


%%

width = 1600;
height = 1000;
x0 = 0;
y0 = 0;

sub_width = 0.21;
sub_height = 0.38;

dy = 0.25;
miny = -3;
maxy = 3;

figure(2)
ax1 = subplot(2,3,4);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-2.e-3,2.e-3,50);
contourf(xax,zax,(Bx_x),levels,'Linecolor','none')
c = colorbar;
c.Ticks = ([-2e-3, -1.e-3, 0.0, 1.0e-3, 2.0e-3]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
caxis([-2.1e-3 2.e-3])
% map = redblue(256);
% colormap(ax1,map(1:256/2,:))
colormap(ax1,redblue)
ylim([-0.03 0.03])
% caxis([-250e1 250e1])
% xlabel('z axis Position (m)')
ylabel('$z$ (cm)','interpreter','latex','Fontsize',30)
xlabel('$x$ (cm)','interpreter','latex','Fontsize',30)
ylabel(c,'$B_x(x,z)$ ($\times 10^{-3}$ T)','interpreter','latex','Fontsize',30)
xlim([0. 0.215])
set(gca, 'XDir','reverse')
set(gca,'Linewidth',1.5)
xticks([0.0, 0.1, 0.2]);
xticklabels([{'0','10','20'}]);
yticks([-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03]);
yticklabels([{'-3','-2','-1','0','1','2','3'}]);
ax1 = get(gca,'Position');
set(gca,'Position',[ax1(1) ax1(2) sub_width sub_height])
% ax = gca();
% c.Ruler.Exponent = 2;
% c.Ruler.TickLabelFormat = '%1.f';


ax2 = subplot(2,3,5);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-2.5e-3,2.5e-3,50);
contourf(xax,zax,(By_x),levels,'Linecolor','none')
colormap(ax2,redblue)
caxis([-2.1e-3 2.e-3])
% xlabel('z axis Position (m)')
xlabel('$x$ (m)','interpreter','latex','Fontsize',30)
ylim([-0.03 0.03])
c = colorbar;
c.Ticks = ([-2e-3, -1.e-3, 0.0, 1.0e-3, 2.0e-3]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
ylabel(c,'$B_y(x,z)$ ($\times 10^{-3}$ T)','interpreter','latex','Fontsize',30)
% ax = gca();
xlim([0. 0.215])
set(gca, 'XDir','reverse')
xticks([0.0, 0.1, 0.2]);
xticklabels([{'0','10','20'}]);
% caxis([-1.1e4 1.1e4])
% c.Ruler.Exponent = 2;
% c.Ruler.TickLabelFormat = '%1.f';
set(gca,'Linewidth',1.5)
set(gca,'ytick',[])
ax2 = get(gca,'Position');
set(gca,'Position',[ax2(1) ax2(2) sub_width sub_height])

ax3 = subplot(2,3,6);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-2.e-3,2.e-3,50);
contourf(xax,zax,(Bz_x),levels,'Linecolor','none')
colormap(ax3,redblue)
xlabel('$x$ (m)','interpreter','latex','Fontsize',30)
% ylabel('x Position (m)')
c = colorbar;
c.Ticks = ([-2e-3, -1.e-3, 0.0, 1.0e-3, 2.0e-3]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
ylim([-0.03 0.03])
ylabel(c,'$B_z(x,z)$ ($\times 10^{-3}$ T)','interpreter','latex','Fontsize',30)
xlim([0. 0.215])
set(gca,'ytick',[])
set(gca, 'XDir','reverse')
xticks([0.0, 0.1, 0.2]);
xticklabels([{'0','10','20'}]);
set(gca,'Linewidth',1.5)
ax3 = get(gca,'Position');
set(gca,'Position',[ax3(1) ax3(2) sub_width sub_height])

ax4 = subplot(2,3,1);
imagesc('XData', actx_BP(:,1), 'YData', acty_BP(1,:), 'CData', plot_temp_bx(:,:,end)')
set(gca, 'XDir','reverse')
colormap(redblue)
ylim([miny-dy/2 maxy+dy/2])
c = colorbar; 
lim = caxis; caxis([-2.0e-3, 2.0e-3])
c.Ticks = ([-2e-3, -1.e-3, 0.0, 1.0e-3, 2.0e-3]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
ylabel(c,'$B_x(x,z)$ ($\times 10^{-3}$ T)','interpreter','latex','Fontsize',30)
set(gca,'xtick',[])
set(gca,'Fontsize',25)
box on
set(gca,'Linewidth',3)
ylabel('$z$ (cm)','interpreter','latex','Fontsize',30)
ax4 = get(gca,'Position');
set(gca,'Position',[ax4(1) ax4(2) sub_width sub_height])

ax5 = subplot(2,3,2);
imagesc('XData', actx_BP(:,1), 'YData', acty_BP(1,:), 'CData', plot_temp_by(:,:,end)')
set(gca, 'XDir','reverse')
colormap(redblue)
ylim([miny-dy/2 maxy+dy/2])
c = colorbar; 
lim = caxis; caxis([-2.0e-3, 2.0e-3])
c.Ticks = ([-2e-3, -1.e-3, 0.0, 1.0e-3, 2.0e-3]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
ylabel(c,'$B_y(x,z)$ ($\times 10^{-3}$ T)','interpreter','latex','Fontsize',30)
set(gca,'Fontsize',25)
set(gca,'ytick',[])
set(gca,'xtick',[])
box on
set(gca,'Linewidth',3)
ax5 = get(gca,'Position');
set(gca,'Position',[ax5(1) ax5(2) sub_width sub_height])

ax6 = subplot(2,3,3);
imagesc('XData', actx_BP(:,1), 'YData', acty_BP(1,:), 'CData', plot_temp_bz(:,:,end)')
set(gca, 'XDir','reverse')
colormap(redblue)
ylim([miny-dy/2 maxy+dy/2])
c = colorbar; 
lim = caxis; caxis([-2.0e-3, 2.0e-3])
c.Ticks = ([-2e-3, -1.e-3, 0.0, 1.0e-3, 2.0e-3]);
c.TickLabels = ([-2.0 -1.0 0.0 1.0 2.0]);
ylabel(c,'$B_z(x,z)$ ($\times 10^{-3}$ T)','interpreter','latex','Fontsize',30)
set(gca,'Fontsize',25)
set(gca,'ytick',[])
set(gca,'xtick',[])
box on
set(gca,'Linewidth',3)
ax6 = get(gca,'Position');
set(gca,'Position',[ax6(1) ax6(2) sub_width sub_height])


%%

xax = zeros(1,npts);
yax = zeros(1,npts);
figure(1)
hold on
xlim([-2.5e-3 2.5e-3])
zlim([-0.0125 0.0125])
ylim([min(zax) max(zax)])
plot3(xax,zax,yax,'k','LineWidth',2)
xlabel('x')
zlabel('y')
ylabel('z position')
box on
% plot_str = [".r",".b",".k"];
% count = 1;
% kk = 1;

for ii=1:npts
    
    figure(1)
    h = quiver3(xax(ii)',zax(ii)',yax(ii)',...
        Ex_x(ii),Ez_x(ii),Ey_x(ii),'k','ShowArrowHead','off','AutoScaleFactor',1.0e-5);
    hs = get(h,'MaxHeadSize');
    set(h,'MaxHeadSize',hs/10)
    plot3(1.0e-5*Ex_x(ii,1), (1.0e-5*Ez_x(ii,1))+zax(1,ii), 1.0e-5*Ey_x(ii,1),'.r')
%     if (ii*dx)<(count*0.3142)
%         plot3(1.0e-5*Ex_x(ii,1), (1.0e-5*Ez_x(ii,1))+zax(1,ii), 1.0e-5*Ey_x(ii,1),plot_str(kk))
%     elseif (ii*dx)>(count*0.3142)
%         count = count + 1;
%         if kk<3
%             kk = kk + 1;
%         elseif kk==3
%             kk = 1;
%         end
%         plot3(1.0e-5*Ex_x(ii,1), (1.0e-5*Ez_x(ii,1))+zax(1,ii), 1.0e-5*Ey_x(ii,1),plot_str(kk))
%     end
     
    pause(0.1)
    delete(h)
    clear h
    hold on

    
end
% quiver3(0,0,zax(npts/4),Ex_x(npts/4),Ey_x(npts/4),Ez_x(npts/4))

%%

x0 = 0;
y0 = 0;
width = 900;
height = 600;

figure(1)
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(0,0.03,50);
contourf(zax,xax,(Emag_x)',levels,'Linecolor','none')
colormap(flipud(gray))
c = colorbar;
xlabel('z Position (m)')
ylabel('x Position (m)')
ylabel(c,'|E|(x,z) (Vm^{-1})','Fontsize',20)
ylim([-0.5 0.5])
c.Ruler.Exponent = -2;
c.Ruler.TickLabelFormat = '%1.1f';

% export_fig('/Volumes/DATA/matlab/wave_verification/real_kx_kz.png',...
%     '-r300')



