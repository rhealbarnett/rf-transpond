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
% define wavenumbers ky and kz (/m); use values given in van eester section IV?? No
% others mentioned
ky = 5.0;
kz = 6.0;

%--
% density N0 (taken from axes on figure 3)
npts = 250;
N0 = logspace(15, 18, npts);

%--
% initialise kx roots arrays, ensure they are complex
kx_arr = zeros(npts, 4);
kx_arr = complex(kx_arr);
kx_roots_arr = zeros(npts, 4);
kx_roots_arr = complex(kx_roots_arr);
kx_coeffs_arr = zeros(npts, 4);
kx_coeffs_arr = complex(kx_coeffs_arr);

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
alpha = 0.5;
beta = 0.5;

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
% k matrix
a11 = ky^2 + kz^2;
a12 = -ky*kx;
a13 = -kz*kx;
a21 = -ky*kx;
a22 = kz^2 + kx^2;
a23 = -ky*kz;
a31 = -kz*kx;
a32 = -ky*kz;
a33 = ky^2 + kx^2;

a = [[a11, a12, a13]
    [a21, a22, a23]
    [a31, a32, a33]];

%-- 
% loop through density values

for ii = 1:npts
    %--
    % electron calcs; density, plasma frequency
    Ne = N0(ii);
    om_pe = sqrt(Ne*e^2/(me*eps0));

    %--
    % ion calcs (95% D, 5% H); density, plasma frequency
    Nd = 0.95*N0(ii);
    om_pd = sqrt(Nd*qd^2/(md*eps0));
    Nh = 0.05*N0(ii);
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
    % rotate cpdt
    cpdt_rot = r'*cpdt*r;
%     cpdt_rot = cpdt;

    %--
    % wave equation rhs
    we_rhs = k0^2*cpdt_rot;

    %--
    % set wave equation to zero to find determinant
    wave_eq = a - we_rhs;

    %--
    % find kx's
    kx_quart = det(wave_eq);
    kx_roots = solve(kx_quart == 0.0, kx);
    kx_roots = double(kx_roots);
    kx_arr(ii,:) = kx_roots;
    
    %--
    % solve using coefficients and roots function
    kx_coeffs = coeffs(kx_quart, 'All');
    kx_coeffs_roots = roots(kx_coeffs);
    kx_roots_arr(ii,:) = kx_coeffs_roots;
    
    %--
    % 'normalise' the coefficients to the coeff on the highest order term
    % put coeffs into the companion matrix, eigenvalues are the roots
    kx_coeffs_norm = kx_coeffs/kx_coeffs(1);

    coeffa = kx_coeffs_norm(2);
    coeffb = kx_coeffs_norm(3);
    coeffc = kx_coeffs_norm(4);
    coeffd = kx_coeffs_norm(5);
    coeffs_matrix = [[0.0, 0.0, 0.0, -coeffd];
                    [1.0, 0.0, 0.0, -coeffc];
                    [0.0, 1.0, 0.0, -coeffb];
                    [0.0, 0.0, 1.0, -coeffa]];
                
    [kx_coeffs_vecs, kx_coeffs_roots] = eig(coeffs_matrix);
    for nn = 1:4
        kx_coeffs_arr(ii,nn) = kx_coeffs_roots(nn,nn);
    end
    
end

%--
% might not be this simple -- kx root values are likely not 'ordered' (MS)
k1 = kx_arr(:,1);
k2 = kx_arr(:,2);
k3 = kx_arr(:,3);
k4 = kx_arr(:,4);

k1_coeff = kx_coeffs_arr(:,1);
k2_coeff = kx_coeffs_arr(:,2);
k3_coeff = kx_coeffs_arr(:,3);
k4_coeff = kx_coeffs_arr(:,4);

k1_root = kx_roots_arr(:,1);
k2_root = kx_roots_arr(:,2);
k3_root = kx_roots_arr(:,3);
k4_root = kx_roots_arr(:,4);

%--
% 'linear' plot axis for N0 powers
N0_plot = linspace(15, 18, npts);

%--
% plot kx's
figure(1)
plot(N0_plot, real(k1_coeff))

hold on

plot(N0_plot, real(k2_coeff))
plot(N0_plot, real(k3_coeff))
plot(N0_plot, real(k4_coeff))
legend('k1', 'k2', 'k3', 'k4')
% ylim([-2.5, 2.5])
xlabel('log10N0 (/m^3)')
ylabel('Real(kx) (/m)')

hold off

figure(2)
plot(N0_plot, imag(k1_coeff))

hold on

plot(N0_plot, imag(k2_coeff))
plot(N0_plot, imag(k3_coeff))
plot(N0_plot, imag(k4_coeff))
legend('k1', 'k2', 'k3', 'k4')
% ylim([-10.0, 10.0])
xlabel('log10N0 (/m^3)')
ylabel('Imag(kx) (/m)')

hold off










