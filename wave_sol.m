%--------------------------------%
% time independent wave solver   %
% V X V X E = k0^2.K.E           %
% rlbarnett c3149416 210917      %
%--------------------------------%

%--
% constants
mu0 = 4.0*pi*1.0e-7;
eps0 = 8.85e-12;
c0 = 1.0/sqrt(eps0*mu0);

%--
% magnetic field (tesla)
B0 = 0.0;

%--
% driver freq
freq = 51.0e6;
om = 2.0*pi*freq;
k0 = om/c0;

%--
% define wavenumbers ky and kz (/m); use values given in van eester section IV?? No
% others mentioned
ky = 5.0;
kz = 6.0;

%--
% density N0 -- set to zero for vacuum case
npts = 10;
xmin = 0.0;
xmax = 2.5;
xax = linspace(xmin,xmax,npts);
dx = (xmax - xmin)/(npts-1);
N0 = 0.0;

%--
% electron constants
e = -1.6022e-19;
me = 9.11e-31;

%--
% ion constants (95% D, 5% H)
qd = abs(e);
mp = 1.67e-27;
md = 2.0*mp;

qh = abs(e);
mh = mp;

%-- 
% cyclotron frequencies
om_ce = e*B0/me;
om_cd = qd*B0/md;
om_ch = qh*B0/mh;

%--
% rotation matrix
alpha = 0.01;
beta = 1.57;

r11 = cos(beta)*cos(alpha);
r12 = cos(beta)*sin(alpha);
r13 = -sin(beta);
r21 = -sin(alpha);
r22 = cos(alpha);
r23 = 0.0;
r31 = sin(beta)*cos(alpha);
r32 = sin(beta)*sin(alpha);
r33 = cos(beta);

r = [[r11, r12, r13]
     [r21, r22, r23]
     [r31, r32, r33]];
 
r_para = r(3,:);
Bvec = B0*r_para;

%--
% electron calcs; density, plasma frequency
Ne = N0;
om_pe = sqrt(Ne*e^2/(me*eps0));

%--
% ion calcs (95% D, 5% H); density, plasma frequency
Nd = 0.95*N0;
om_pd = sqrt(Nd*qd^2/(md*eps0));
Nh = 0.05*N0;
om_ph = sqrt(Nh*qh^2/(mh*eps0));

%--
% cold plasma dielectric tensor elements
s = 1.0 - om_pe^2/(om^2 - om_ce^2) - om_pd^2/(om^2 - om_cd^2) - om_ph^2/(om^2 - om_ch^2);
d = om_ce*om_pe^2/(om*(om^2 - om_ce^2)) + om_cd*om_pd^2/(om*(om^2 - om_cd^2)) + om_ch*om_ph^2/(om*(om^2 - om_ch^2));
p = 1.0 - om_pe^2/om^2 - om_pd^2/om^2 - om_ph^2/om^2;

%--
% cold plasma delectric tensor
cpdt = [[s, -1i*d, 0.0]
        [1i*d, s, 0.0]
        [0.0, 0.0, p]];
    
%--
% set up rhs matrix (multiples E field)

waveeq_mat = zeros(3.0*npts, 9.0*npts);

for ii=2,3.0*npts:3;
    eq1_row = 3.0*ii - 2.0;
    eq2_row = 3.0*ii - 1.0;
    eq3_row = 3.0*ii;
    for kk=1,npts;
        
    
    waveeq_mat(eq1_row,ii) = -1i*ky;
    waveeq_mat(eq1_row,ii-1) = 0.0;
    waveeq_mat(eq1_row,ii+1) = -1i*kz;
    waveeq_mat(eq2_row,ii) = -1.0;
    waveeq_mat(eq2_row,ii-1) = -1i*ky*(dx/2.0);
    waveeq_mat(eq2_row,ii+1) = 0.0;
    waveeq_mat(eq3_row,ii) = 0.0;
    waveeq_mat(eq3_row,ii-1) = -1i*kz*(dx/2.0);
    waveeq_mat(eq3_row,ii+1) = -1.0;
    
    

























