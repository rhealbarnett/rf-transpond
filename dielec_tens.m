%-----------------------------------------%
% Populate cold plasma dielectric tensor  %
% arbitrary species                       %
% rlbarnett c3149416, 121217              %
%-----------------------------------------%

%------------------------------------------------------------------------
% input for mass, density can be vectors
% output will be plasma frequency for each mass 

function [om_c,om_p,cpdt,s_arr,d_arr,p_arr] = dielec_tens(q,B0,n,m,om,eps0,npts)

    cpdt = zeros(3,3,npts);
    s_arr = zeros(1,npts);
    d_arr = zeros(1,npts);
    p_arr = zeros(1,npts);
    
    msize = size(m);

    for ii=1:msize(1)
        om_p(ii,:) = plasma_freq(q,n,m(ii,:),eps0);
        om_c(ii,:) = cyclo_freq(q,B0,m(ii,:));
    end

    s = 1.0 - sum((om_p.^2)./(om^2 - om_c.^2),1);
    d = sum(om_c.*om_p.^2./(om*(om^2 - om_c.^2)),1);
    p  = 1.0 - sum((om_p.^2/om^2),1);

    s_arr(1,:) = s;
    d_arr(1,:) = d;
    p_arr(1,:) = p;

    %--
    % cold plasma delectric tensor
    cpdt(1,1,:) = s;
    cpdt(2,2,:) = s;
    cpdt(1,2,:) = -1i*d;
    cpdt(2,1,:) = 1i*d;
    cpdt(3,3,:) = p;

%     cpdt_arr(:,:,nn) = rot'*cpdt_arr(:,:,nn)*rot;

end