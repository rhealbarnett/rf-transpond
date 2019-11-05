%-----------------------------------------%
% Populate cold plasma dielectric tensor  %
% arbitrary species                       %
% rlbarnett c3149416, 121217              %
%-----------------------------------------%

%------------------------------------------------------------------------
% input for mass, density can be vectors
% output will be plasma frequency for each mass 

function [om_c,om_p,cpdt,s_arr,d_arr,p_arr,sigma] = dielec_tens(q,B0,n,m,om,eps0,npts)

    cpdt = zeros(3,3,npts);
    s_arr = zeros(1,npts);
    d_arr = zeros(1,npts);
    p_arr = zeros(1,npts);
    
    msize = size(m);
    nsize = size(n);
    
%     dampFac = 1.0;
%     np_bound = floor(0.2*npts);
%     ax = linspace(0,pi,np_bound);
%     damp0 = (cos(ax)+1)/2;
%     damp = ones(1,npts);
%     damp(1:np_bound) = damp(1:np_bound) + dampFac*1i*damp0;
%     damp(end-np_bound+1:end) = damp(end-np_bound+1:end) + dampFac*1i*fliplr(damp0);
    om = om*ones(1,npts);
%     om = om.*damp;

    for ii=1:msize(1)
        if nsize(1)~=1
            om_p(ii,:) = plasma_freq(q(ii,:),n(ii,:),(m(ii,:)),eps0);
        else
            om_p(ii,:) = plasma_freq(q(ii,:),n,(m(ii,:)),eps0);
        end
        om_c(ii,:) = cyclo_freq(q(ii,:),B0,(m(ii,:)));
    end

    s = 1.0 - sum(((om_p).^2)./((om).^2 - (om_c).^2),1);
    d = sum((om_c).*(om_p).^2./((om).*((om).^2 - (om_c).^2)),1);
    p  = 1.0 - sum(((om_p).^2./(om).^2),1);

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
    
    sigma = zeros(size(cpdt));
    
    for ii=1:npts
        sigma(:,:,ii) = 1i*om(1,ii)*(eps0*eye(3) - cpdt(:,:,ii));
    end

%     cpdt_arr(:,:,nn) = rot'*cpdt_arr(:,:,nn)*rot;

end