xdata = h5read('../../lapd_calcs/LAPD_Rhea/ReEx.hdf','/page1.axes1.image1/xdata');
ydata = h5read('../../lapd_calcs/LAPD_Rhea/ReEx.hdf','/page1.axes1.image1/ydata');

Rex_data = h5read('../../lapd_calcs/LAPD_Rhea/ReEx.hdf','/page1.axes1.image1/cdata');
Iex_data = h5read('../../lapd_calcs/LAPD_Rhea/ImEx.hdf','/page1.axes1.image1/cdata');

exdata = Rex_data + Iex_data*1i;

ex_mag = abs(exdata);

figure(1)
levels = linspace(-30,30,100);
contourf(xdata,ydata,ex_mag',levels,'Linecolor','none')
set(gca,'colorscale','log')
colorbar
lim = caxis;
caxis([min(levels), max(levels)]);
cbh=colorbar('v');
set(cbh,'YTick',[-30:3:30])
ylabel(cbh,'E_x (units??)')


slice0_index = find(ydata>=0);
slice0_index = slice0_index(1);
yval = ydata(slice0_index);

slice1_index = find(ydata>=0.1);
slice1_index = slice1_index(1);
yval1 = ydata(slice1_index);

slice2_index = find(ydata>=0.2);
slice2_index = slice2_index(1);
yval2 = ydata(slice2_index);

slice3_index = find(ydata>=0.3);
slice3_index = slice3_index(1);
yval3 = ydata(slice3_index);

slice4_index = find(ydata>=0.4);
slice4_index = slice4_index(1);
yval4 = ydata(slice4_index);

datamax = max([ex_mag(:,slice0_index);ex_mag(:,slice1_index);...
    ex_mag(:,slice2_index);ex_mag(:,slice3_index);ex_mag(:,slice4_index)]);
datamin = min([ex_mag(:,slice0_index);ex_mag(:,slice1_index);...
    ex_mag(:,slice2_index);ex_mag(:,slice3_index);ex_mag(:,slice4_index)]);

figure(2)
subplot(5,1,5)
plot(xdata,ex_mag(:,slice0_index),'-o')
ylabel('E_x(y=0)','Fontsize',16)
ax = gca;
ax.YAxis.Exponent = 1;
% ylim([datamin, datamax])

subplot(5,1,4)
plot(xdata,ex_mag(:,slice1_index))
ylabel('E_x(y=0.1)','Fontsize', 16)
ax = gca;
ax.YAxis.Exponent = -2;
set(gca,'xtick',[])
% ylim([datamin, datamax])

subplot(5,1,3)
plot(xdata,ex_mag(:,slice2_index))
ylabel('E_x(y=0.2)','Fontsize', 16)
ax = gca;
ax.YAxis.Exponent = -3;
set(gca,'xtick',[])
% ylim([datamin, datamax])

subplot(5,1,2)
plot(xdata,ex_mag(:,slice3_index))
ylabel('E_x(y=0.3)','Fontsize', 16)
ax = gca;
ax.YAxis.Exponent = -3;
set(gca,'xtick',[])
% ylim([datamin, datamax])

subplot(5,1,1)
plot(xdata,ex_mag(:,slice4_index))
ylabel('E_x(y=0.4)','Fontsize', 16)
ax = gca;
ax.YAxis.Exponent = 1;
set(gca,'xtick',[])
% ylim([datamin, datamax])

