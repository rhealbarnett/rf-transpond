%-------------------------------------------------------------------------%
% Projection plots. 
% rlbarnett c3149416 200109
%-------------------------------------------------------------------------%

width = 600;
height = 900;
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
xlabel('z axis Position (m)')
ylabel('x axis Position (m)')
ylabel(c,'E_x(x,z) (Vm^{-1})','Fontsize',20)

ax2 = subplot(3,1,2);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-1.0e3,1.0e3,50);
contourf(zax,xax,real(Ey_x)',levels,'Linecolor','none')
colormap(ax2,redblue)
xlabel('z axis Position (m)')
ylabel('x axis Position (m)')
c = colorbar;
ylabel(c,'E_y(x,z) (Vm^{-1})','Fontsize',20)

ax3 = subplot(3,1,3);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-1.8,1.8,50);
contourf(zax,xax,real(Ez_x)',levels,'Linecolor','none')
colormap(ax3,redblue)
xlabel('z axis Position (m)')
ylabel('x axis Position (m)')
c = colorbar;
ylabel(c,'E_z(x,z) (Vm^{-1})','Fontsize',20)

export_fig('/Volumes/DATA/LAPD/matlab/wave_projection/ez_xz.png',...
    '-r300')