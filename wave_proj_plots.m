%-------------------------------------------------------------------------%
% Projection plots. 
% rlbarnett c3149416 200109
%-------------------------------------------------------------------------%

width = 600;
height = 600;
x0 = 0;
y0 = 0;

figure(1)
ax1 = subplot(3,1,1);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-215,215,50);
contourf(zax,xax,real(Ex_x)',levels,'Linecolor','none')
c = colorbar;
% map = redblue(256);
% colormap(ax1,map(1:256/2,:))
colormap(ax1,redblue)
caxis([-215 215])
% xlabel('z axis Position (m)')
ylabel('x Position (m)')
ylabel(c,'E_x(x,z) (Vm^{-1})','Fontsize',20)
ylim([-0.4 0.401])
% ax = gca();
c.Ruler.Exponent = 2;
c.Ruler.TickLabelFormat = '%1.f';
set(gca,'xtick',[])

ax2 = subplot(3,1,2);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-1.e3,1.e3,50);
contourf(zax,xax,real(Ey_x)',levels,'Linecolor','none')
colormap(ax2,redblue)
% xlabel('z axis Position (m)')
ylabel('x Position (m)')
c = colorbar;
ylabel(c,'E_y(x,z) (Vm^{-1})','Fontsize',20)
% ax = gca();
ylim([-0.4 0.401])
caxis([-1.e3 1.e3])
c.Ruler.Exponent = 2;
c.Ruler.TickLabelFormat = '%1.f';
set(gca,'xtick',[])

ax3 = subplot(3,1,3);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-2.0,2.0,50);
contourf(zax,xax,real(Ez_x)',levels,'Linecolor','none')
colormap(ax3,redblue)
xlabel('z Position (m)')
ylabel('x Position (m)')
c = colorbar;
ylabel(c,'E_z(x,z) (Vm^{-1})','Fontsize',20)
ylim([-0.4 0.401])
% c.Ruler.Exponent = 3;
caxis([-2.0 2.0])
c.Ruler.TickLabelFormat = '%1.f';

export_fig('/Volumes/DATA/LAPD/matlab/wave_projection/e_xz_reducedxaxis.png',...
    '-r300')