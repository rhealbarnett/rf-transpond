%------------------------------------------%
% poisson solver                           %
% adapted from python code 150617          %
% rlbarnett c3149416 191017                %
%------------------------------------------%

%function static_pot = poisson_sol(Ne, Ni, lamby, lambz, e, eps0, dx, npts)

%--
% initialise coefficient matrix and right hand side vector
coeff_mat = zeros(npts,3.0*npts);
rhs = zeros(npts,1);

%--
% calculate charge density
rho = (dx^2*e/eps0)*(N0e - N0i);

%--
% counter
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
    
    %--
    % rolling boundary conditions
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
    
    %--
    % populate coefficient matrix
    coeff_mat(ii,iixm) = 1.0;
    coeff_mat(ii,iiym) = 0.0;
    coeff_mat(ii,iizm) = 0.0;
    coeff_mat(ii,iix) = -2.0;
    coeff_mat(ii,iiy) = dx^2*lamby^2;
    coeff_mat(ii,iiz) = dx^2*lambz^2;
    coeff_mat(ii,iixp) = 1.0;
    coeff_mat(ii,iiyp) = 0.0;
    coeff_mat(ii,iizp) = 0.0;
    
    %--
    % populate right hand side vector
    rhs(ii) = rho(ii);
    
    kk = kk + 3;
    
end

%--
% set simple PEC boundaries for now
rhs(1) = 0.0;
rhs(npts) = 0.0;

%--
% reduce matrix size
coeff_mat = sparse(coeff_mat);

%--
% solve Ax = b
static_pot = coeff_mat\rhs;

static_potx = static_pot(1:3:3*npts);
static_poty = static_pot(2:3:3*npts);
static_potz = static_pot(3:3:3*npts);

%end
    
    







