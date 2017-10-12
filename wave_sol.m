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
freq = 51.0e7;
om = 2.0*pi*freq;
nu = 0.05*om;
om = om + (1i*nu);
k0 = om/c0;
wavel0 = (2*pi)/k0;

%--
% define wavenumbers ky and kz (/m); use values given in van eester section IV?? No
% others mentioned
ky = 0.0;
kz = 0.0;

%--
% spatial domain
npts = 1000;
dx = 0.025;
xmin = 0.0;
% xmax = 10.0;
xmax = npts*dx;
% npts = ((xmax - xmin)/dx);
xax = linspace(xmin,xmax,npts);
% dx = (xmax - xmin)/(npts - 1);

%--
% density -- set to zero for vacuum case
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
    
% cpdt_rot = r'*cpdt*r;
    
%--
% set up rhs matrix (multiples E field)

% syms ky kz k0 dx
% cpdt = sym('R%d%d',[3,3]);

waveeq_mat = zeros(3*npts, 3*npts);
% waveeq_mat = sym('O%d%d', [3*npts,3*npts]);
ii = 1;

for eq1=1:3:3*npts
    
    eq2 = eq1+1;
    eq3 = eq2+1;
    
    iiex = ii;
    iiey = ii+1;
    iiez = ii+2;
    iiexm = iiex - 3;
    iieym = iiey - 3;
    iiezm = iiez - 3;
    iiexp = iiex + 3;
    iieyp = iiey + 3;
    iiezp = iiez + 3;

    %--
    % set up periodic boundary conditions
    if ((iiexm) & (iieym) & (iiezm)) <= 0
        iiexm = 3*npts + iiexm;
        iieym = 3*npts + iieym;
        iiezm = 3*npts + iiezm;
    end
    
    %--
    % this loop doesn't work if it's set up like the previous one???
    if ((iiexp)) >= (3*npts)
        iiexp = iiexp - 3*npts;
        iieyp = iieyp - 3*npts;
        iiezp = iiezp - 3*npts;
    end          
    
    %--
    % fill matrix
    waveeq_mat(eq1,iiexm) = 0.0;
    waveeq_mat(eq1,iieym) = -1i*ky;
    waveeq_mat(eq1,iiezm) = -1i*kz;
    waveeq_mat(eq1,iiex) = 2.0*dx*(ky^2 + kz^2 - k0^2*cpdt(1,1));
    waveeq_mat(eq1,iiey) = -2.0*dx*k0^2*cpdt(1,2);
    waveeq_mat(eq1,iiez) = -2.0*dx*k0^2*cpdt(1,3);
    waveeq_mat(eq1,iiexp) = 0.0;
    waveeq_mat(eq1,iieyp) = 1i*ky;
    waveeq_mat(eq1,iiezp) = 1i*kz;
    
    waveeq_mat(eq2,iiexm) = -1i*ky*(dx/2.0);
    waveeq_mat(eq2,iieym) = -1.0;
    waveeq_mat(eq2,iiezm) = 0.0;
    waveeq_mat(eq2,iiex) = -dx^2*k0^2*cpdt(2,1);
    waveeq_mat(eq2,iiey) = dx^2*(kz^2 - k0^2*cpdt(2,2)) + 2.0;
    waveeq_mat(eq2,iiez) = -dx^2*(ky*kz + k0^2*cpdt(2,3));
    waveeq_mat(eq2,iiexp) = 1i*ky*(dx/2.0);
    waveeq_mat(eq2,iieyp) = -1.0;
    waveeq_mat(eq2,iiezp) = 0.0;
    
    waveeq_mat(eq3,iiexm) = -1i*kz*(dx/2.0);
    waveeq_mat(eq3,iieym) = 0.0;
    waveeq_mat(eq3,iiezm) = -1.0;
    waveeq_mat(eq3,iiex) = -dx^2*k0^2*cpdt(3,1);
    waveeq_mat(eq3,iiey) = -dx^2*(ky*kz + k0^2*cpdt(3,2));
    waveeq_mat(eq3,iiez) = dx^2*(ky^2 - k0^2*cpdt(3,3)) + 2.0;
    waveeq_mat(eq3,iiexp) = 1i*kz*(dx/2.0);
    waveeq_mat(eq3,iieyp) = 0.0;
    waveeq_mat(eq3,iiezp) = -1.0;
    
    ii = ii + 3;

end

%--
% set up rhs vector
rhs = zeros(3*npts,1);
rhs((3*npts/2)+1) = 0.0;
rhs((3*npts/2)+2) = 0.0;
rhs((3*npts/2)+3) = 1.0;

%--
% calculation solution as waveeq_mat^-1*rhs
%-- COMMENT OUT IF DOING SYMBOLIC MATRIX --%
sol = (waveeq_mat)\rhs;

%----------------------plots-----------------------%
figure(1)
plot(xax, real(sol(1:3:3*npts)))

hold on

plot(xax, real(sol(2:3:3*npts)))
plot(xax, real(sol(3:3:3*npts)))
xlabel('Position')
ylabel('Amplitude (?)')
legend('Re[Ex] (?)', 'Re[Ey] (?)', 'Re[Ez] (?)')

hold off

figure(2)
plot(xax, imag(sol(1:3:3*npts)))

hold on

plot(xax, imag(sol(2:3:3*npts)))
plot(xax, imag(sol(3:3:3*npts)))
xlabel('Position')
ylabel('Amplitude (?)')
legend('Im[Ex] (?)', 'Im[Ey] (?)', 'Im[Ez] (?)')

hold off


























