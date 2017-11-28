%------------------------------------------%
% calculate perturbed velocity             %
% eq A2 DVE 2015                           %
% rlbarnett c3149416 281117                %
%------------------------------------------%

denom = 1.0 ./ (om^2 - om_ce.^2);
cross_term = zeros(npts,3);
dot_term = zeros(npts,3);

for ii=1:npts
    
    cross_term(ii,:) = cross(ee_vec(ii,:),Bvec);
    dot_term(ii,:) = (1i / om)*producte(ii)*om_ce_vec;
    
end
    

v1e = denom*(-1i*om.*ee_vec + cross_term + dot_term);