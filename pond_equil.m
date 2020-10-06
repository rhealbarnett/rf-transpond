%-------------------------------------------------------------------------%
% pond_equil : predict final density profile for ponderomotive force
%              given initial density and steady state electric field.
%
% rlbarnett, 2020-10-06
%
% 
%
%
%

function n_final = pond_equil(E,n_init,omega,v_thermal)

    shape = omega./E;
    n_final = n_init.*exp(-(1./2.)*(shape./v_thermal).^2);

end