%------------------------------------------%
% poisson solver                           %
% adapted from python code 150617          %
% rlbarnett c3149416 191017                %
%------------------------------------------%

%function static_pot = poisson_sol(Ne, Nh, Nd, lamby, lambz, e, eps0, dx, npts)

coeff_mat = zeros(npts,3.0*npts);
rhs = zeros(npts,1);

rho = (dx^2*e/eps0)*(Ne - (Nh + Nd));

kk = 1;

for ii=1:npts
    
    iix = kk;
    iiy = kk + 1;
    iiz = kk + 2;
    iixm = iix - 3;
    iiym = iiy - 3;
    iizm = iiz - 3;
    iixp = iix + 3;
    iiyp = iiy + 3;
    iizp = iiz + 3;
    
    
    if ((iixp)) > 3.0*npts
        iixp = iixp - 3.0*npts;
        iiyp = iiyp - 3.0*npts;
        iizp = iizp - 3.0*npts;
    end
    if ((iixm) & (iiym) & (iizm)) <= 0
        iixm = 3.0*npts + iixm;
        iiym = 3.0*npts + iiym;
        iizm = 3.0*npts + iizm;
    end
    
    coeff_mat(ii,iixm) = 1.0;
    coeff_mat(ii,iiym) = 0.0;
    coeff_mat(ii,iizm) = 0.0;
    coeff_mat(ii,iix) = -2.0;
    coeff_mat(ii,iiy) = dx^2*lamby^2;
    coeff_mat(ii,iiz) = dx^2*lambz^2;
    coeff_mat(ii,iixp) = 1.0;
    coeff_mat(ii,iiyp) = 0.0;
    coeff_mat(ii,iizp) = 0.0;
    
    rhs(ii) = rho;
    
    kk = kk + 3;
    
end

coeff_mat = sparse(coeff_mat);

static_pot = coeff_mat\rhs;

%end
    
    







