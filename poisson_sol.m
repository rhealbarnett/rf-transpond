%------------------------------------------%
% poisson solver                           %
% adapted from python code 150617          %
% rlbarnett c3149416 191017                %
%------------------------------------------%

%--
% constants
mu0 = 4.0*pi*1.0e-7;
eps0 = 8.85e-12;
c0 = 1.0/sqrt(eps0*mu0);
e = -1.6022e-19;
qd = abs(e);
qh = abs(e);

%--
% spatial domain
npts = 10;
dx = 0.02;
xmin = 0.0;
xmax = 0.2;
xax = linspace(xmin, xmax, npts);

%--
% density
N0 = 5.0e17;
Ne = N0;
Nd = 0.95*N0;
Nh = 0.15*N0;
Ni = Nd + Nh;

%--
% charge density
rho = -(e/eps0)*(Ne - Ni);

%--
% initialise potential array and set up simple boundary conditions
phi = zeros(npts);
phi(1) = 0.0;
phi(end) = 0.0;
phi_old = zeros(npts);

%--
% solver parameters
threshold = 1.0e-5;
nmax = 10000;

%--
% solve

for ii=nmax
    
    phi(2:end-1) = (1.0/2.0)*(phi(1:end-2) - phi(3:end) - dx^2*rho);
    
    rms_err = sqrt(sum((phi - phi_old)^2)/npts);
    
    phi_old(:) = phi;
    
    if rms_err<=threshold
        break
    else
        continue
    end
end










