% ----------------------plots----------------------- %
figure(16)
% set(gca,'YTickLabel','%.2f')

subplot(4,1,1)
plot(nxax,om*const.mu0*source,'r')
ylabel('Source (i\omega\mu_0J_{y,z})')
xlim([xmin,xmax])

subplot(4,1,2)
plot(nxax, real(rf_ex), 'k')
ylabel('E_x (Vm^{-1})')

hold on

plot(nxax, imag(rf_ex), '--')
set(gca, 'XTickLabel', [])
legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
xlim([xmin,xmax])

hold off

subplot(4,1,3)
plot(nxax, real(rf_ey), 'k')
ylabel('E_y (Vm^{-1})')

hold on

plot(nxax, imag(rf_ey), '--')
set(gca, 'XTickLabel', [])
legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
xlim([xmin,xmax]);

hold off

subplot(4,1,4)
plot(nxax, real(rf_ez), 'k')
ylabel('E_z (Vm^{-1})')

hold on

plot(nxax, imag(rf_ez), '--')
xlabel('Position (m)')
legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
xlim([xmin,xmax])

hold off




