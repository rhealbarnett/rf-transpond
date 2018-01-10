%------------------------------------------%
% calculate ponderomotive acceleration     %
% expression from DVE 2015                 %
% rlbarnett c3149416 131117                %
%------------------------------------------%

rf_E = zeros(3,npts);
rf_E(1,:) = rf_ex;
rf_E(2,:) = rf_ey;
rf_E(3,:) = rf_ez;

om_ce_vec = zeros(3,npts);
om_ch_vec = zeros(3,npts);
om_cd_vec = zeros(3,npts);
ee_vec = zeros(3,npts);
eh_vec = zeros(3,npts);
ed_vec = zeros(3,npts);

for ii=1:npts
    
    om_ce_vec(:,ii) = e*Bvec(:,ii) / me(ii);
    om_ch_vec(:,ii) = qh*Bvec(:,ii) / mh(ii);
    om_cd_vec(:,ii) = qd*Bvec(:,ii) / md(ii);
    
    ee_vec(:,ii) = e*rf_E(:,ii) / me(ii);
    eh_vec(:,ii) = qh*rf_E(:,ii) / mh(ii);
    ed_vec(:,ii) = qd*rf_E(:,ii) / md(ii);

end

ee_mag = sqrt(ee_vec(1,:).*conj(ee_vec(1,:)) + ee_vec(2,:).*conj(ee_vec(2,:)) + ee_vec(3,:).*conj(ee_vec(3,:)));
eh_mag = sqrt(eh_vec(1,:).*conj(eh_vec(1,:)) + eh_vec(2,:).*conj(eh_vec(2,:)) + eh_vec(3,:).*conj(eh_vec(3,:)));
ed_mag = sqrt(ed_vec(1,:).*conj(ed_vec(1,:)) + ed_vec(2,:).*conj(ed_vec(2,:)) + ed_vec(3,:).*conj(ed_vec(3,:)));

producte = dot(ee_vec,om_ce_vec);
productd = dot(ed_vec,om_cd_vec);
producth = dot(eh_vec,om_ch_vec);

% producte = zeros(3,npts);
% producth = zeros(3,npts);
% productd = zeros(3,npts);

% for ii=1:npts
%     
%     producte(:,ii) = dot(ee_vec(:,ii),om_ce_vec);
%     producth(:,ii) = dot(eh_vec(:,ii),om_ch_vec);
%     productd(:,ii) = dot(ed_vec(:,ii),om_cd_vec(:,ii);
%     
% end

pond_pote = (1.0 / 4.0)*(1.0 ./ (om^2 - om_ce.^2)).*(ee_mag.^2 - abs(producte / om).^2 - (2.0*om_ce / om).*imag(conj(ee_vec(1,:).*ee_vec(2,:))));
pond_poth = (1.0 / 4.0)*(1.0 ./ (om^2 - om_ch.^2)).*(eh_mag.^2 - abs(producth / om).^2 - (2.0*om_ch / om).*imag(conj(eh_vec(1,:).*eh_vec(2,:))));
pond_potd = (1.0 / 4.0)*(1.0 ./ (om^2 - om_cd.^2)).*(ed_mag.^2 - abs(productd / om).^2 - (2.0*om_cd / om).*imag(conj(ed_vec(1,:).*ed_vec(2,:))));
pond_poti = pond_poth + pond_potd;

a_pondex = -gradient(pond_pote,dx);
a_pondix = -gradient(pond_poti,dx);
a_pondhx = -gradient(pond_poth,dx);
a_ponddx = -gradient(pond_potd,dx);

a_pondey = -lamby*pond_pote;
a_pondez = -lambz*pond_pote;

a_pondiy = -lamby*pond_poti;
a_pondiz = -lambz*pond_poti;
a_pondhy = -lamby*pond_poth;
a_pondhz = -lambz*pond_poth;
a_ponddy = -lamby*pond_potd;
a_ponddz = -lambz*pond_potd;



