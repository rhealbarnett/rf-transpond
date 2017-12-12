%-----------------------------------------%
% Populate cold plasma dielectric tensor  %
% arbitrary species                       %
% rlbarnett c3149416, 121217              %
%-----------------------------------------%


function ans = om_c(npts,q,B0,m)

    ans = q*B0./m;
    
end

function ans = om_p(npts,q,N,eps0,m)

    ans = sqrt(q^2*N./(eps0*m));
    
end

function ans = cpdt(npts,N,B0,om,q,m,sign,eps0)

    om_c = q*B0./m;
    om_p = sqrt(q^2*N./(eps0*m));
    
end
    