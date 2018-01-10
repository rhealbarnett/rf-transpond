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

accele = zeros(3,npts);
acceli = zeros(3,npts);

accele(1,:) = accelex;
accele(2,:) = acceley;
accele(3,:) = accelez;

acceli(1,:) = accelix;
acceli(2,:) = acceliy;
acceli(3,:) = acceliz;

v_drifte = zeros(3,npts);
v_drifti = zeros(3,npts);

om_cimag = sqrt(om_ch.^2 + om_cd.^2);

for ii=1:npts
    
    v_drifte(:,ii) = (1.0/om_ce(ii)).*cross(accele(:,ii),Bvec(:,ii));
    v_drifti(:,ii) = (1.0/om_cimag(ii)).*cross(acceli(:,ii),Bvec(:,ii));
    
end

vd_perp1e = v_drifte(1,:);
vd_perp2e = v_drifte(2,:);

vd_perp1i = v_drifti(1,:);
vd_perp2i = v_drifti(2,:);
