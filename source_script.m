% ------------------------------------------------------------------- %
% call file for different E field values in wave solver / transport 
% 250219 rlbarnett, c3149416
% ------------------------------------------------------------------- %

npts = 2048;
density_pert = sparse(npts,10);
n_arr = sparse(npts,10);
Evec_init = zeros(1,10);
iter = 1;


for source_mult = 1000:4000:37000
    
    source_mult
    
    transport_1d;
    Efin = max(abs(sqrt(Efield)));
    Evec_init(1,iter) = Emag;
    Evec_fin(1,iter) = Efin;
    pert = n_init - n_new;
    density_pert(:,iter) = interp1(nxax,pert,xax,'linear');
    n_arr(:,iter) = n_new;
    
    iter = iter + 1;
    
    
end

%%

levels = linspace(-20e15,20e15,50);
% levels = logspace(0,log(3.0e16),50);

c = colormap('redblue');

dp11 = sign(density_pert).*log10(abs(density_pert));

figure(42)
contourf(xax,Evec_fin,density_pert',levels,'Linecolor','none')
ylabel('Max RF |E_{||}| (Vm^{-1})','Fontsize',16)
xlabel('Position (m)','Fontsize',16)
title('Density perturbation (m^{-3})')
colormap(gca,c);
% set(gca,'colorscale','linear')
lim = caxis;
caxis([-20e15 20e15])
colorbar;
    