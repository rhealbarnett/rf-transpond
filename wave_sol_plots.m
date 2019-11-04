% ----------------------plots----------------------- %
figure(16)
% set(gca,'YTickLabel','%.2f')

subplot(4,1,1)
plot(xax,ey_source,'b','Linewidth',2)
ylabel('Source (i\omega\mu_0J_{y})')
set(gca, 'XTickLabel', [])
xlim([xmin,xmax])
set(gca,'Fontsize',20)

subplot(4,1,2)
plot(xax, real(rf_ex), 'k','Linewidth',2)
ylabel('E_x (Vm^{-1})')

hold on

plot(xax, imag(rf_ex), '--r','Linewidth',2)
set(gca, 'XTickLabel', [])
legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
xlim([xmin,xmax])
set(gca,'Fontsize',20)

hold off

subplot(4,1,3)
plot(xax, real(rf_ey), 'k','Linewidth',2)
ylabel('E_y (Vm^{-1})')

hold on

plot(xax, imag(rf_ey), '--r','Linewidth',2)
set(gca, 'XTickLabel', [])
legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
xlim([xmin,xmax]);
set(gca,'Fontsize',20)

hold off

subplot(4,1,4)
plot(xax, real(rf_ez), 'k','Linewidth',2)
ylabel('E_z (Vm^{-1})')

hold on

plot(xax, imag(rf_ez), '--r','Linewidth',2)
xlabel('Position (m)')
legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
xlim([xmin,xmax])
set(gca,'Fontsize',20)

hold off




