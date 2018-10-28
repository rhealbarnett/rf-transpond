%-----------------------------------------%
% Populate cold plasma dielectric tensor  %
% arbitrary species                       %
% rlbarnett c3149416, 121217              %
%-----------------------------------------%

function [om_c,om_p,cpdt] = dielec_tens(q,B0,n,m,om,eps0,npts)

cpdt = zeros(3,3,npts);
s_arr = zeros(1,npts);
d_arr = zeros(1,npts);
p_arr = zeros(1,npts);
    
for nn=1:npts
    
    for ii=1:length(m)
        om_p(ii,:) = plasma_freq(q,n,m(ii),eps0);
        om_c(ii,1) = cyclo_freq(q,B0,m(ii));
    end
    
    s = 1.0 - sum((om_p.^2)./(om^2 - om_c.^2));
    d = sum((om_c*om_p.^2/(om*om^2 - om_c.^2)));
    p  = 1.0 - sum(om_p.^2/om^2);
    
%     s = 1.0 - om_pe(nn)^2/(om^2 - om_ce(nn)^2) - om_pd(nn)^2/(om^2 - om_cd(nn)^2) - om_ph(nn)^2/(om^2 - om_ch(nn)^2);
%     d = om_ce(nn)*om_pe(nn)^2/(om*(om^2 - om_ce(nn)^2)) + om_cd(nn)*om_pd(nn)^2/(om*(om^2 - om_cd(nn)^2)) +...
%         om_ch(nn)*om_ph(nn)^2/(om*(om^2 - om_ch(nn)^2));
%     p = 1.0 - om_pe(nn)^2/om^2 - om_pd(nn)^2/om^2 - om_ph(nn)^2/om^2;

    s_arr(1,nn) = s;
    d_arr(1,nn) = d;
    p_arr(1,nn) = p;
    
    %--
    % cold plasma delectric tensor
    cpdt(1,1,nn) = s;
    cpdt(2,2,nn) = s;
    cpdt(1,2,nn) = -1i*d;
    cpdt(2,1,nn) = 1i*d;
    cpdt(3,3,nn) = p;
    
%     cpdt_arr(:,:,nn) = rot'*cpdt_arr(:,:,nn)*rot;
   
end

r_arr = s_arr + d_arr;
l_arr = s_arr - d_arr;

end