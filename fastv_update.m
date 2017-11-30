%------------------------------------------%
% calculate perturbed velocity             %
% eq A2 DVE 2015                           %
% rlbarnett c3149416 281117                %
%------------------------------------------%

denome = 1.0 ./ (om^2 - om_ce.^2);
denomi = 1.0 ./ (om^2 - om_cimag.^2);
om_ci_vec = sqrt(om_cd_vec.^2 + om_ch_vec.^2);
ei_vec = sqrt(eh_vec.^2 + ed_vec.^2);
producti = sqrt(producth.^2 + productd.^2);
cross_terme = zeros(3,npts);
dot_terme = zeros(3,npts);
cross_termi = zeros(3,npts);
dot_termi = zeros(3,npts);

for ii=1:npts
    
    cross_terme(:,ii) = cross(ee_vec(:,ii),Bvec);
    dot_terme(:,ii) = (1i / om)*producte(ii)*om_ce_vec;
    cross_termi(:,ii) = cross(ei_vec(:,ii),Bvec);
    dot_termi(:,ii) = (1i / om)*producti(ii)*om_ci_vec;
    
end
    

v1e = denome*(-1i*om.*ee_vec + cross_terme + dot_terme);
v1i = denomi*(-1i*om.*ei_vec + cross_termi + dot_termi);