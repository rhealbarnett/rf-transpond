%-------------------------------------------------------------------------%
% Projection plots. 
% rlbarnett c3149416 200109
%-------------------------------------------------------------------------%

width = 900;
height = 600;
x0 = 0;
y0 = 0;

figure(1)
set(gcf,'Position',[x0 y0 width height],'Color','w')
ax1 = figure(1);
levels = linspace(-1.8,1.8,50);
contourf(zax,xax,real(Ez_x)',levels,'Linecolor','none')
colormap(ax1,redblue)
c = colorbar;