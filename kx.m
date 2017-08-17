%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% values/parameters from van eester 2015  %
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

clear all; close all

syms kx; 

%--
% constants
mu0 = 4.0*pi*1.0e-7;
eps0 = 8.85e-12;
c0 = 1.0/sqrt(eps0*mu0);

%--
% magnetic field (tesla)
B0 = 2.5;

%--
% driver freq
freq = 51.0e6;
om = 2.0*pi*freq;
k0 = om/c0;

%--
% define wavenumbers ky and kz (/m); use values given in section IV?? No
% others mentioned
ky = 5.0;
kz = 6.0;

%--
% density N0 (taken from axes on figure 3)
N0 = logspace(15, 18, 1000);
count = length(N0);

%--
% initialise kx roots array
kx_arr = zeros(count, 4);

%--
% electron constants
e = -1.6022e-19;
me = 9.11e-31;

%--
% ion constants (assume 100% D plasma for now)
qd = 2.0*abs(e);
mp = 1.67e-27;
md = 2.0*mp;

%-- 
% cyclotron frequencies
om_ce = e*B0/me;
om_cd = qd*B0/md;

%-- 
% loop through density values

for ii = 1:count
    syms z;
    %--
    % electron calcs
    Ne = N0(ii);
    om_pe = sqrt(Ne*e^2/(me*eps0));

    %--
    % ion calcs
    Nd = N0(ii);
    om_pd = sqrt(Nd*qd^2/(md*eps0));

    %--
    % cold plasma dielectric tensor elements
    s = 1.0 - om_pe.^2/(om^2 - om_ce^2) - om_pd.^2/(om^2 - om_cd^2);
    d = om_ce*om_pe.^2/(om*(om^2 - om_ce^2)) + om_cd*om_pd.^2/(om*(om^2 - om_cd^2));
    p = 1.0 - om_pe.^2/om^2 - om_pd.^2/om^2;

    %--
    % cold plasma delectric tensor
    cpdt = [[s, -j*d, 0.0]
            [j*d, s, 0.0]
            [0.0, 0.0, p]];

    %--
    % rotation matrix
    alpha = 0.01;
    beta = 0.01;

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

    %-- 
    % rotate cpdt
    cpdt_rot = transpose(r)*cpdt*r;

    %--
    % wave equation rhs
    we_rhs = k0^2*cpdt_rot;

    %-- 
    % k matrix
    a11 = ky^2 + kz^2;
    a12 = - ky*kx;
    a13 = -kz*kx;
    a21 = -ky*kx;
    a22 = kz^2 + kx^2;
    a23 = -ky*kz;
    a31 = -kz*kx;
    a32 = -ky*kz;
    a33 = ky^2 + kx^2;

    % %-- 
    % % k - k0**2.K matrix
    % a11 = ky^2 + kz^2 - k0^2*s;
    % a12 = k0^2*j*d - ky*kx;
    % a13 = -kz*kx;
    % a21 = -ky*kx - k0^2*j*d;
    % a22 = kz^2 + kx^2 - k0^2*s;
    % a23 = -ky*kz;
    % a31 = -kz*kx;
    % a32 = -ky*kz;
    % a33 = ky^2 + kx^2 - k0^2*p;

    a = [[a11, a12, a13]
        [a21, a22, a23]
        [a31, a32, a33]];

    %--
    % set wave equation to zero to find determinant
    wave_eq = a - we_rhs;

    %--
    % find kx's
    kx_quart = det(wave_eq);
    kx_roots = solve(kx_quart == 0, kx);
    kx_roots = eval(kx_roots);
    
    %--
    % getting solution in terms of 'free' variable z. More eqns than
    % unknowns? Should be solutions for any value of z-- set to 1?
    z = 1;  

    kx_arr(ii,:) = kx_roots;
    
end

k1 = kx_arr(:,1);
k2 = kx_arr(:,2);
k3 = kx_arr(:,3);
k4 = kx_arr(:,4);


semilogx(N0, real(k1))

hold on

semilogx(N0, real(k2))
semilogx(N0, real(k3))
semilogx(N0, real(k4))
legend('k1', 'k2', 'k3', 'k4')
xlim([min(N0), max(N0)])
xlabel('log10N0 (/m^3)')
ylabel('Real(kx) (/m)')

hold off

















