%------------------------------------------------------------------------%
% ion_visc : calculate ion collision time and para viscosity coefficient. 
% 
% rlbarnett 20210513
%
% INPUTS
% n : ion density (m^-3)
% Ti : ion temperature (eV)
% m : ion mass (kg)
% Z : atomic number
% 
% OUTPUTS
% tau : ion collision time
% eta : parallel ion viscosity coefficient
%
%------------------------------------------------------------------------%

function [tau, eta] = ion_visc(n, Ti, m, Z)

    const = constants();

    % Log of the Coulomb algorithm (usually 10-20)
    coulomb_log = 17;
    
    tau = (((4.0*pi*const.eps0)^2)*3.0*sqrt(m)*(Ti*const.e)^(3/2))./...
        (4.0*sqrt(pi)*coulomb_log*const.e^4*Z^4*n);
    
    eta = 0.96*(Ti*const.e)*n.*tau;

end

