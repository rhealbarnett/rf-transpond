%----------------------------------------------------------------------%
% calculate poynting flux from wave solver
% rlbarnett c3149416
% 20/06/07
%----------------------------------------------------------------------%

%----------------------------------------------------------------------%
% INPUTS
% rf_e(x,y,z): rf e field solution from wave_sol.m calculation.
% k(x,y): wavenumbers used as inputs in wave_sol.m function to calculate
%         rf e solution. 
% ax: spatial axis.
% om: rf driving frequency. 
%
% OUPUTS
% poyn: the poynting flux Sx, Sy and Sz, in Wm^-2. 
%----------------------------------------------------------------------%

function poyn = poynting(rf_ex, rf_ey, rf_ez, kx, ky, ax, om)

    const = constants();
    mu0 = const.mu0;
    npts = length(ax);

    [rf_bx, rf_by, rf_bz] = B_RF(ax,kx,ky,om,rf_ex,rf_ey,rf_ez,0);
    
    rf_e = [rf_ex; rf_ey; rf_ez];
    rf_h = (1.0/mu0)*[rf_bx; rf_by; rf_bz];
    
    poyn = zeros(3,npts);
    
    for ii=1:npts
    
        poyn(:,ii) = 0.5*real(cross(rf_e(:,ii),conj(rf_h(:,ii))));
        
    end

end