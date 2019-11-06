% ----------------------------------------------------------------------- %
% Script to read and plot matlab (.mat) files for the coupled 1D 
% wave / transport solver. 
%
% rlbarnett c3149416 191106
% ----------------------------------------------------------------------- %

filepath = '/Volumes/DATA/LAPD/matlab/results_jsource_kyzero/';

for ii=106
    
    filename = strcat(filepath, 'coupled_transport_', num2str(ii),'.mat');

    load(filename)
    
    figure(1)
%     plot(nxax,n_init,'color',[0,0,0]+0.8,'Linewidth',4)
    hold on
    plot(nxax,n_new,'k','Linewidth',2)
    xlabel('Position (m)')
    ylabel('Density (m^{-3})')
    set(gca,'Fontsize',20)
    xlim([min(nxax) max(nxax)])
    
    figure(2)
%     plot(vxax,vx_init,'color',[0,0,0]+0.8,'Linewidth',4)
    hold on
    plot(vxax,vx_new,'k','Linewidth',2)
    
    figure(3)
    plot(xax,abs(rf_ez).^2)
    hold on
    
    figure(4)
    plot(xax,rf_ez)
    hold on
    
end
