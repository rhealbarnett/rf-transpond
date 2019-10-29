% ----------------------------------------------------------------------- %
% Function to calculate Ey driven by Jy antenna source.
% rlbarnett c3149416 191030.
%
% Uses near field antenna calculation (see Cheng eq 11-18b):
% 
%           E_theta = [p/(4*pi*eps0*R^3)]*sin(theta).
%
% At the midplane y = 0, theta = 0. R is calculated using the distance
% from the antenna to the location of the solution domain and the length
% of the domain. 
%
% Inputs
% xax: solution domain in metres. 
% dipole_mom: the dipole moment p. Used to scale the electric field output
% to the desired value. 
% source_dist: distance from the current source to the domain in metres.
%
% Outputs
% e_source: electric field calculated at each location in the domain. 
% ----------------------------------------------------------------------- %


function [etheta_source] = e_source(xax,dipole_mom,source_dist)

    const = constants();
    eps0 = const.eps0;
    
    R = sqrt(source_dist^2 + xax.^2);
    etheta_source = dipole_mom./(4.0*pi*eps0.*(R.^3));
    
end
    