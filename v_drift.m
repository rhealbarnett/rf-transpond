%------------------------------------------%
% calculate perpendicular drift velocity   %
% v_drift = (1.0/OMc)*(a X Bvec)           %
% rlbarnett c3149416 161117                %
%------------------------------------------%

accelex = pressex + static_ex + a_pondex;
acceley = pressey + static_ey + a_pondey;
accelez = pressez + static_ez + a_pondez;

accelix = pressix + static_ex + a_pondix;
acceliy = pressiy + static_ey + a_pondiy;
acceliz = pressiz + static_ez + a_pondiz;

accele = zeros(npts,3);
acceli = zeros(npts,3);

accele(:,1) = accelex;
accele(:,2) = acceley;
accele(:,3) = accelez;

acceli(:,1) = accelix;
acceli(:,2) = acceliy;
acceli(:,3) = acceliz;

v_drifte = zeros(npts,3);

for ii=1:npts
    
    v_drifte(ii,:) = (1.0/om_ce)*cross(accele(ii,:),Bvec);
    
end

vd_perp1e = v_drifte(:,1);
vd_perp2e = v_drifte(:,2);
