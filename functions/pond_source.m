%------------------------------------------------------------------%
% Function calculates ponderomotive acceleration
% rlbarnett c3149416 updated 191205.
%
% INPUTS
% component: string input, options are 'para', 'perp' or 'total'. This 
%            will set which component of the electric field is used in the 
%            calculation of the parallel ponderomotive acceleration. If 
%            desired output is a single array of spatial values, i.e. 
%            summed over species and components, include second input in
%            cell array e.g. {'total',1}
% Efield: {ex, ey, ez}. Can enter empty string if a particular component
%         of the field is not being used e.g. {ex, ey, ''}. Must be on 
%         uniform grid. 
% m: mass. If multiple species, enter as column array e.g. [me; mi] for
%    electron and ion mass. 
% q: charge. If multiple species, enter as column array e.g. [e; q] for
%    electron and ion charge. Position in the array needs to be consistent
%    with position in the mass array.
% om_cyc: cyclotron frequency. If multiple species, enter as column array 
%         e.g. [om_ce; om_ci] for electron and ion frequencies. Position 
%         in the array needs to be consistent with position in charge and 
%         mass arrays.
% omega: driving frequency for RF field. 
% npts: number of points in the spatial domain. 
% dz: grid spacing. Needs to be an array of size (1,npts). If dz is 
%     constant over the domain, enter as same value at each location.
%
% OUPUTS
% ans: array of ponderomotive acceleration for input into an eq of 
%      motion, the velocity equation. 2x2xnpts array of values. 
%      Parallel component is stored in row index 1, 
%      perpendicular in row index 2. Column index is dictated by 
%      index of species input. E.g. for m = [me; mi], q = [e; q] and 
%      om_cyc = [om_ce; om_ci], 'total' output is
%
%       2x2xnpts array
%       ans(1,1,:) = PA Epara,e    ans(1,2,:) = Epara,i
%       ans(2,1,:) = PA Epara,e    ans(2,2,:) = Epara,i
%
% 
%
%------------------------------------------------------------------%

function [Ediff,pf] = pond_source(component,Efield,m,q,om_cyc,omega,dz,mix,damping)

    rf_ex = Efield{1};
    rf_ey = Efield{2};
    rf_ez = Efield{3};
    
    npts = length(rf_ex);
    ax = linspace(-4,4,npts);
    
    if mix
        Emix = imag((conj(rf_ey).*...
                    rf_ex)-(conj(rf_ex).*rf_ey));
    else
        Emix = zeros(1,length(rf_ex));
    end
    
    Esquared = [abs(Efield{1}).^2; abs(Efield{2}).^2; abs(Efield{3}).^2; Emix];
    Esize = size(Esquared);
    Ediff = zeros(numel(Efield)+1,length(rf_ex));
       
%     for ii = 2:length(dz)
        for jj = 1:Esize(1)
%             Ediff(jj,ii-1) = -(dz(1,ii)/(dz(1,ii-1)*(dz(1,ii) + dz(1,ii-1))))*Esquared(jj,ii-1) +...
%                 ((dz(1,ii) - dz(1,ii-1))/(dz(1,ii)*dz(1,ii-1)))*Esquared(jj,ii) +...
%                 (dz(1,ii-1)/(dz(1,ii)*(dz(1,ii) + dz(1,ii-1))))*Esquared(jj,ii+1);
            Ediff(jj,2:end-1) = (Esquared(jj,3:end) - Esquared(jj,1:end-2)) ./ (2.0*dz);
            Ediff(jj,1) = (-3*Esquared(jj,1) + 4*Esquared(jj,2) - Esquared(jj,3))/(2.0*dz);
            Ediff(jj,end) = (3*Esquared(jj,end) - 4*Esquared(end-1) + Esquared(end-2))/(2.0*dz);
        end       
%     end
    
    if strcmp(component{1},'total')
        
        para = 1;
        perp = 1;
        
    end

    msize = size(m);
    
    pf_const = zeros(msize(1),1);
    pond = zeros(2,2,length(rf_ex));
    
    for ii=1:msize(1)
    
        pf_const(ii,1) = q(ii,1)^2/(4.0*m(ii,1));
        
        
        if strcmp(component{1},'para')

            pond(1,ii,:) = (pf_const(ii,1)/omega^2).*Ediff(3,:);
                
        elseif strcmp(component{1},'perp')

            pond(2,ii,:) = (pf_const(ii,1)).*(((Ediff(1,:) + Ediff(2,:))/(omega^2 - om_cyc(ii,1)^2))...
                    + (om_cyc(ii,1)/(omega*(omega^2 - om_cyc(ii,1)^2))).*Ediff(4,:));
                
        elseif strcmp(component{1},'total')
            
             pond(1,ii,:) = (pf_const(ii,1)/omega^2).*Ediff(3,:);
             pond(2,ii,:) = (pf_const(ii,1)).*(((Ediff(1,:) + Ediff(2,:))/(omega^2 - om_cyc(ii,1)^2))...
                    + (om_cyc(ii,1)/(omega*(omega^2 - om_cyc(ii,1)^2))).*Ediff(4,:));

        end
        
    end
    
    if damping{1}
        dampSize = floor(damping{2}*npts);
        dampLoc = ax(dampSize);
        ind = find(damping{3}==dampLoc);
        for ii=1:msize(1)
%             ind = find(damping{2}<=-2.4);
            for jj=1:msize(1)
                pond(ii,jj,1:ind) = 0.0;
                pond(ii,jj,length(damping{3})-ind-1:end) = 0.0;
            end
            
        end
    end
        
    
    pf = (1.0/m(2,1)).*pond;
    
    if component{2}
        intermediate = sum(pf,1);
        pf = squeeze(sum(intermediate,2))';
    else
    end
%     ans = (1.0/m(2,1)).*pond;
%     ans = pond;

end