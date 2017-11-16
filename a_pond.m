%------------------------------------------%
% calculate ponderomotive acceleration     %
% expression from DVE 2015                 %
% rlbarnett c3149416 131117                %
%------------------------------------------%

rf_E = zeros(npts,3);
rf_E(:,1) = rf_ex;
rf_E(:,2) = rf_ey;
rf_E(:,3) = rf_ez;

om_ce_vec = e*Bvec / me;
om_ch_vec = qh*Bvec / mh;
om_cd_vec = qd*Bvec / md;

ee_vec = e*rf_E / me;
eh_vec = qh*rf_E / mh;
ed_vec = qd*rf_E / md;

ee_mag = sqrt(ee_vec(:,1).*conj(ee_vec(:,1)) + ee_vec(:,2).*conj(ee_vec(:,2)) + ee_vec(:,3).*conj(ee_vec(:,3)));
eh_mag = sqrt(eh_vec(:,1).*conj(eh_vec(:,1)) + eh_vec(:,2).*conj(eh_vec(:,2)) + eh_vec(:,3).*conj(eh_vec(:,3)));
ed_mag = sqrt(ed_vec(:,1).*conj(ed_vec(:,1)) + ed_vec(:,2).*conj(ed_vec(:,2)) + ed_vec(:,3).*conj(ed_vec(:,3)));

% om_ce_mag = sqrt(om_ce_vec(:,1).*conj(om_ce_vec(:,1)) + om_ce_vec(:,2).*conj(om_ce_vec(:,2)) + om_ce_vec(:,3).*conj(om_ce_vec(:,3)));
% om_ch_mag = sqrt(om_ch_vec(:,1).*conj(om_ch_vec(:,1)) + om_ch_vec(:,2).*conj(om_ch_vec(:,2)) + om_ch_vec(:,3).*conj(om_ch_vec(:,3)));
% om_cd_mag = sqrt(om_cd_vec(:,1).*conj(om_cd_vec(:,1)) + om_cd_vec(:,2).*conj(om_cd_vec(:,2)) + om_cd_vec(:,3).*conj(om_cd_vec(:,3)));

producte = zeros(npts,1);
producth = zeros(npts,1);
productd = zeros(npts,1);

for ii=1:npts
    
    producte(ii) = dot(ee_vec(ii,:),om_ce_vec);
    producth(ii) = dot(eh_vec(ii,:),om_ch_vec);
    productd(ii) = dot(ed_vec(ii,:),om_cd_vec);
    
end

pond_pote = (1.0 / 4.0)*(1.0 / (om^2 - om_ce^2))*(ee_mag.^2 - abs(producte / om).^2 - (2.0*om_ce / om)*imag(conj(ee_vec(:,1).*ee_vec(:,2))));
pond_poth = (1.0 / 4.0)*(1.0 / (om^2 - om_ch^2))*(eh_mag.^2 - abs(producth / om).^2 - (2.0*om_ch / om)*imag(conj(eh_vec(:,1).*eh_vec(:,2))));
pond_potd = (1.0 / 4.0)*(1.0 / (om^2 - om_cd^2))*(ed_mag.^2 - abs(productd / om).^2 - (2.0*om_cd / om)*imag(conj(ed_vec(:,1).*ed_vec(:,2))));
pond_poti = pond_poth + pond_potd;

a_pondex = gradient(pond_pote,dx);
a_pondix = gradient(pond_poti,dx);

a_pondey = -lamby*pond_pote;
a_pondez = -lambz*pond_pote;

a_pondiy = -lamby*pond_poti;
a_pondiz = -lambz*pond_potz;



