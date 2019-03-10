% ------------------------------------------------------------------- %
% call file for different E field values in wave solver / transport 
% 250219 rlbarnett, c3149416
% ------------------------------------------------------------------- %

npts = 2048;
density_pert = sparse(npts,10);
Evec = zeros(1,10);
iter = 1;


for source_mult = 1000:2000:25000
    
    source_mult
    
    transport_1d;
    Evec(1,iter) = Emag;
    pert = n_init - n_new;
    density_pert(:,iter) = interp1(nxax,pert,xax,'linear');
    
    iter = iter + 1;
    
    
end

%%

levels = linspace(-30e15,30e15,20);

c = colormap('redblue');

figure(42)
contourf(xax,Evec,density_pert',levels,'Linecolor','none')
ylabel('Max RF Re[Ex] (Vm^{-1})')
xlabel('Position (m)')
title('Density perturbation (m^{-3})')
colormap(gca,c);
lim = caxis;
caxis([-30e15 30e15])
colorbar;
    