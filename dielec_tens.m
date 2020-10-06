%-----------------------------------------%
% Populate cold plasma dielectric tensor  %
% arbitrary species                       %
% rlbarnett c3149416, 121217              %
%-----------------------------------------%

%------------------------------------------------------------------------%
% INPUTS
% q: charge, as a column vector for each species e.g. q = [-e; e..] for 
%    however many species in the plasma, in C. 
% B0: background magnetic field, in T, as a scalar.
% n: density in m^-3 over spatial domain, column for each species (must be
%    ordered the same as q). 
% m: species mass in kg, column for each species, ordered as with q and n.
% om: RF driver frequency in Hz, as a scalar.
% eps0: vacuum permittivity, in SI units.
% npts: number of spatial grid points.
% damping: switch 0 or 1 to turn real(sigma) damping profile off or on.
%
% OUTPUTS:
% om_c: cyclotron frequency.
% om_p: plasma frequency.
% cpdt: cold plasma dielectric tensor, 3X3 matrix at each spatial location.
% s_arr, d_arr, p_arr: cold plasma dielectric elements S, P and D at each      
%                      spatial location.
% sig: sigma tensor, 3X3 at each spatial location.
%------------------------------------------------------------------------%


%------------------------------------------------------------------------
% input for mass, density can be vectors
% output will be plasma frequency for each mass 

function [om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q,B0,n,m,om,eps0,npts,damping)

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
%     om = om*ones(1,npts) + 0.1i*om*ones(1,npts);
%     m = m + 0.1i*m;

    for ii=1:msize(1)
        if nsize(1)~=1
            om_p(ii,:) = plasma_freq(q(ii,:),n(ii,:),(m(ii,:)),eps0);
        else
            om_p(ii,:) = plasma_freq(q(ii,:),n,(m(ii,:)),eps0);
        end
        om_c(ii,:) = cyclo_freq(q(ii,:),B0,(m(ii,:)));
    end

    s = 1.0 - sum((om.*(om_p).^2)./(real(om).*(om.^2 - (om_c).^2)),1);
    d = sum((om_c).*(om_p).^2./(real(om).*((om).^2 - (om_c).^2)),1);
    p = 1.0 - sum(((om_p).^2./(real(om).*om)),1);

    
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
    
    sig = zeros(size(cpdt));
    
    for ii=1:npts
        sig(:,:,ii) = 1i*om(1,ii)*(eps0*eye(3) - cpdt(:,:,ii));
    end

%     cpdt_arr(:,:,nn) = rot'*cpdt_arr(:,:,nn)*rot;

    if damping{1}
        
        const = constants();
        eps0 = const.eps0;

        np_bound = floor(damping{2}*npts);
        ax = linspace(0,pi,np_bound);
        damp0 = (cos(ax)+1)/2;
        damp_sig = ones(1,npts);
        damp_sig(1:np_bound) = damp0*-1 + max(damp0);
        damp_sig(end-np_bound+1:end) = fliplr(damp_sig(1:np_bound));
        shift11 = max(imag(sig(1,1,:)));
        shift33 = max(imag(sig(3,3,:)));
        
        for ii=1:npts
            
            sig(1,1,ii) = sig(1,1,ii)*damp_sig(1,ii);
            sig(1,1,ii) = sig(1,1,ii) + ((damp_sig(1,ii)*-1 + 1)*shift11)*1.e3;
%             sig(1,1,ii) = sig(1,1,ii) + (1i*sig(1,1,ii)*damp_sig(1,ii));
%             sig(2,2,ii) = damp_sig(1,ii)*(sig(2,2,ii));
%             sig(2,2,ii) = sig(2,2,ii) + (1i*(sig(2,2,ii)*(damp_sig(1,ii)))+shift11)*1.0e3;
            sig(2,2,ii) = sig(1,1,ii);
            sig(3,3,ii) = damp_sig(1,ii)*(sig(3,3,ii));
            sig(3,3,ii) = sig(3,3,ii) + ((damp_sig(1,ii)*-1 + 1)*shift33)*0.5e1;
            cpdt(1,1,ii) = eps0 - (1.0/(1i*om(1,ii)))*sig(1,1,ii);
            cpdt(2,2,ii) = cpdt(1,1,ii);%eps0 - (1.0/(1i*om(1,ii)))*sig(2,2cd ph,ii);
            cpdt(3,3,ii) = eps0 - (1.0/(1i*om(1,ii)))*sig(3,3,ii);
            
        end
        
    end

end







