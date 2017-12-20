%-----------------------------------------%
% Populate cold plasma dielectric tensor  %
% arbitrary species                       %
% rlbarnett c3149416, 121217              %
%-----------------------------------------%

function ans = cpdt(npts,N,B0,om,q,m,sign,eps0)

    om_c = q*B0./m;
    om_p = sqrt(q^2*N./(eps0*m));
    
end
    