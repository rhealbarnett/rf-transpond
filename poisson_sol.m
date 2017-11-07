%------------------------------------------%
% poisson solver                           %
% adapted from python code 150617          %
% rlbarnett c3149416 191017                %
%------------------------------------------%

coeff_mat = zeros(npts,npts);
rhs = zeros(npts,1);

mult = dx^2*(lamby^2 + lambz^2) - 2.0;
rho = (dx^2*e/eps0)*(Ne - (Nh + Nd));

for ii=1:npts
    
    iip = ii + 1;
    iim = ii - 1;
    
    if iip > npts
        iip = iip - npts;
    end
    if iim <= 0
        iim = npts + iim;
    end
    
    coeff_mat(ii,ii) = mult;
    coeff_mat(ii,iim) = 1.0;
    coeff_mat(ii,iip) = 1.0;
    
    rhs(ii) = rho(ii);
    
end

sol = coeff_mat\rhs;
    
    







