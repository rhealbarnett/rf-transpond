%------------------------------------------------------------------%
% Function calculates ponderomotive acceleration
% rlbarnett c3149416 updated 191127.
%
% INPUTS
% m: mass of the particle being accelerated.
% mp: proton mass from the velocity update equation.
% omega: driving frequency for the fast E field.
% q: charge of the particle being accelerated.
% Efield: absolute value of the fast electric field squared, the quantity
%       to be differentiated.
% dx: array of dx values. 
% varargin: option to include cyclotron frequency required for calculation
%       of the ion ponderomotive acceleration.
%
% OUPUTS
% ans: the ponderomotive acceleration for input into an eq of motion, the 
%       velocity equation.
%
%------------------------------------------------------------------%

function [ans] = pond_source(component,Efield,m,q,om_cyc,omega,npts,dz)

    Esquared = [Efield{1}.^2; Efield{2}.^2; Efield{3}.^2];
    Esize = size(Esquared);
    Ediff = zeros(numel(Efield),npts);
    
    for jj = 1:Esize(1)
        for ii = 2:length(dz)
            Ediff(jj,ii-1) = -(dz(1,ii)/(dz(1,ii-1)*(dz(1,ii) + dz(1,ii-1))))*Esquared(jj,ii-1) +...
                ((dz(1,ii) - dz(1,ii-1))/(dz(1,ii)*dz(1,ii-1)))*Esquared(jj,ii) +...
                (dz(1,ii-1)/(dz(1,ii)*(dz(1,ii) + dz(1,ii-1))))*Esquared(jj,ii+1);
        end
    end
    
    if strcmp(component,'full')
        
        para = 1;
        perp = 1;
        
    end

    msize = size(m);
    
    pf_const = zeros(msize(1),1);
    pond = {[],msize(1)};
%     pond_para = zeros(msize,npts);
%     pond_perp = zeros(msize,npts);
    
    for ii=1:msize(1)
    
        pf_const(ii,1) = q(ii,1)^2/(4.0*m(ii,1));
        
        
        if strcmp(component,'para') || para

                pond{1,ii} = (pf_const(ii,1)/omega^2).*Ediff(3,:);
                
        elseif strcmp(component,'perp') || perp

                pond{2,ii} = (pf_const(ii,1)).*(((Ediff(1,:) + Ediff(2,:))/(omega^2 - om_cyc(1,1)^2))...
                    + (om_cyc(ii,1)/(omega*(omega^2 - om_cyc(ii,1)^2)))*imag((conj(Efield{2})*Efield{1})...
                    -(conj(Efield{1})*Efield{2})));

        end
        
    end
        

    ans = (1.0/m(2,1))*pond;
end