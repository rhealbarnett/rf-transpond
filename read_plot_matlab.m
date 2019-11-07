% ----------------------------------------------------------------------- %
% Script to read and plot matlab (.mat) files for the coupled 1D 
% wave / transport solver. 
%
% rlbarnett c3149416 191106
% ----------------------------------------------------------------------- %

equib = load('/Volumes/DATA/LAPD/matlab/lapd_equib_superrefined.mat','vx_new','n_new','xmin',...
    'xmax','npts','cs');
filepath = '/Volumes/DATA/LAPD/matlab/results_jsource_kyzero/';

mach_init = equib.vx_new/equib.cs;
n_init = equib.n_new;
xax = linspace(equib.xmin,equib.xmax,equib.npts);

clear equib

for ii=22154
    
    filename = strcat(filepath, 'coupled_transport_', num2str(ii),'.mat');

    load(filename)
    
    figure(1)
    plot(nxax,n_init,'color',[0,0,0]+0.8,'Linewidth',4)
    hold on
    box on
    plot(nxax,n_new,'k','Linewidth',2)
    xlabel('Position (m)')
    ylabel('Density (m^{-3})')
    xlim([min(nxax) max(nxax)])
    
    figure(2)
    plot(vxax,mach_init,'color',[0,0,0]+0.8,'Linewidth',4)
    hold on
    box on
    plot(vxax,vx_new/cs,'k','Linewidth',2)
    xlabel('Position (m)')
    ylabel('Mach #')
    xlim([min(vxax) max(vxax)])
    
    
    figure(3)
    plot(xax,rf_ez,'k','Linewidth',2)
    hold on
    xlabel('Position (m)')
    ylabel('RF E_z (Vm^{-1})')
    xlim([min(xax) max(xax)])
    
    figure(4)
    plot(xax,abs(rf_ez).^2,'k','Linewidth',2)
    hold on
    xlabel('Position (m)')
    ylabel('|RF E_z (Vm^{-1})|^2')
    xlim([min(xax) max(xax)])

    
end
