
function [density_perturbation, denpert_avg] = den_pert(n_init, n_final, oscil)
    
    density_perturbation = (n_init - n_final) ./ n_init;
    
    if oscil
        
        [yupper,ylower] = envelope(density_perturbation);
        denpert_avg = yupper + ylower;
        
    end
    
   

end