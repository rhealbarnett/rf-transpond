% ----------------------plots----------------------- %
figure(16)
% set(gca,'YTickLabel','%.2f')

subplot(4,1,1)
plot(xax,om*const.mu0*source,'b')
ylabel('Source (i\omega\mu_0J_{y,z})')
xlim([xmin,xmax])

subplot(4,1,2)
plot(xax, real(rf_ex), 'k')
ylabel('E_x (Vm^{-1})')

hold on

plot(xax, imag(rf_ex), '--r')
set(gca, 'XTickLabel', [])
legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
xlim([xmin,xmax])

hold off

subplot(4,1,3)
plot(xax, real(rf_ey), 'k')
ylabel('E_y (Vm^{-1})')

hold on

plot(xax, imag(rf_ey), '--r')
set(gca, 'XTickLabel', [])
legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
xlim([xmin,xmax]);

hold off

subplot(4,1,4)
plot(xax, real(rf_ez), 'k')
ylabel('E_z (Vm^{-1})')

hold on

plot(xax, imag(rf_ez), '--r')
xlabel('Position (m)')
legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
xlim([xmin,xmax])

hold off




