%------------------------------------------------------------------%
% Function calculates ponderomotive acceleration
% rlbarnett c3149416 updated 191205.
%
% INPUTS
% component: string input, options are 'para', 'perp' or 'total'. This 
%            will set which component of the electric field is used in the 
%            calculation of the parallel ponderomotive acceleration.
% Efield: {ex, ey, ez}. Can enter empty string if a particular component
%         of the field is not being used e.g. {ex, ey, ''}.
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
% ans: cell array of ponderomotive acceleration for input into an eq of 
%      motion, the velocity equation. 2x2 cell, array of spatial values in
%      each cell. Parallel component is stored in row index 1, 
%      perpendicular in row index 2. Returns 0x0 double array in cell if 
%      particular component is not calculated. Column index is dictated by 
%      index of species input. E.g. for m = [me; mi], q = [e; q] and 
%      om_cyc = [om_ce; om_ci], 'total' output is
%
%       2x2 cell array
%       {PA Epara, e 1xnpts double}  {PA Epara, i 1xnpts double}
%       {PA Eperp, e 1xnpts double}  {PA Eperp, i 1xnpts double}
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
    
    if strcmp(component,'total')
        
        para = 1;
        perp = 1;
        
    end

    msize = size(m);
    
    pf_const = zeros(msize(1),1);
    pond = {};
    
    for ii=1:msize(1)
    
        pf_const(ii,1) = q(ii,1)^2/(4.0*m(ii,1));
        
        
        if strcmp(component,'para')

            pond{1,ii} = (pf_const(ii,1)/omega^2).*Ediff(3,:);
                
        elseif strcmp(component,'perp')

            pond{2,ii} = (pf_const(ii,1)).*(((Ediff(1,:) + Ediff(2,:))/(omega^2 - om_cyc(1,1)^2))...
                    + (om_cyc(ii,1)/(omega*(omega^2 - om_cyc(ii,1)^2)))*imag((conj(Efield{2}).*Efield{1})...
                    -(conj(Efield{1}).*Efield{2})));
                
        elseif strcmp(component,'total')
            
             pond{1,ii} = (pf_const(ii,1)/omega^2).*Ediff(3,:);
             pond{2,ii} = (pf_const(ii,1)).*(((Ediff(1,:) + Ediff(2,:))/(omega^2 - om_cyc(1,1)^2))...
                    + (om_cyc(ii,1)./(omega*(omega^2 - om_cyc(ii,1)^2))).*imag((conj(Efield{2}).*Efield{1})...
                    -(conj(Efield{1}).*Efield{2})));

        end
        
    end
        

%     ans = (1.0/m(2,1))*pond;
    ans = pond;
end