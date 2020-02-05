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
levels = linspace(-251,251,50);
contourf(zax,yax,(Ex_y)',levels,'Linecolor','none')
c = colorbar;
% map = redblue(256);
% colormap(ax1,map(1:256/2,:))
colormap(ax1,redblue)
caxis([-250 250])
% xlabel('z axis Position (m)')
ylabel('y Position (m)')
ylabel(c,'E_x(y,z) (Vm^{-1})','Fontsize',20)
ylim([-0.4 0.401])
% ax = gca();
c.Ruler.Exponent = 2;
c.Ruler.TickLabelFormat = '%1.f';
set(gca,'xtick',[])

ax2 = subplot(3,1,2);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-1.1e3,1.1e3,50);
contourf(zax,yax,(Ey_y)',levels,'Linecolor','none')
colormap(ax2,redblue)
% xlabel('z axis Position (m)')
ylabel('y Position (m)')
c = colorbar;
ylabel(c,'E_y(y,z) (Vm^{-1})','Fontsize',20)
% ax = gca();
ylim([-0.4 0.401])
caxis([-1.1e3 1.1e3])
c.Ruler.Exponent = 2;
c.Ruler.TickLabelFormat = '%1.f';
set(gca,'xtick',[])

ax3 = subplot(3,1,3);
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = linspace(-2.02,2.02,50);
contourf(zax,yax,(Ez_y)',levels,'Linecolor','none')
colormap(ax3,redblue)
xlabel('z Position (m)')
ylabel('y Position (m)')
c = colorbar;
ylabel(c,'E_z(y,z) (Vm^{-1})','Fontsize',20)
ylim([-0.4 0.401])
% c.Ruler.Exponent = 3;
caxis([-2.0 2.0])
c.Ruler.TickLabelFormat = '%1.f';

% export_fig('/Volumes/DATA/LAPD/matlab/wave_projection/e_xz_reducedxaxis.png',...
%     '-r300')


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



