
function [density_perturbation,mean_pert] = den_pert(n_init, n_final, nxax, plots)
    
    density_perturbation = (n_init - n_final) ./ n_init;
   
    ax = find(nxax >= -0.5 & nxax <= 0.5);
    wind_np = length(ax);
    uni_ax = linspace(min(nxax(ax)),max(nxax(ax)),wind_np);
    
    wind = hann(wind_np);
    pert_wind = interp1(uni_ax,wind,nxax(ax),'linear');
    
    density_wind = pert_wind.*density_perturbation(ax);
    
    if plots
        
        thesis_fig(1,nxax(ax),density_perturbation(ax),'$R_n$','Position (m)',...
            1.5,[0 0 0]+0.7,'')
        hold on
        thesis_fig(1,nxax(ax),density_wind,'$R_n$','Position (m)',...
            1.5,'k','')
        xlim([-0.5 0.5])
        plot(-0.06*ones(1,wind_np),linspace(min(density_perturbation(ax)),...
            max(density_perturbation(ax)),wind_np),'--b')
        plot(0.06*ones(1,wind_np),linspace(min(density_perturbation(ax)),...
            max(density_perturbation(ax)),wind_np),'--b')
        
    end
    
    mean_pert = mean(density_wind);
   
end