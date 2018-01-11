%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% values/parameters from van eester 2015  %
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

syms kx

%--
% initialise kx roots arrays, ensure they are complex
kx_arr = zeros(npts, 4);
kx_arr = complex(kx_arr);

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
    % wave equation rhs
    we_rhs = k0^2*cpdt;

    %--
    % set wave equation to zero to find determinant
    wave_eq = a - we_rhs(:,:,ii);

    %--
    % find kx's
    kx_quart = det(wave_eq);
    kx_coeffs = coeffs(kx_quart, 'All');
    kx_coeffs_roots = roots(kx_coeffs);
    kx_arr(ii,:) = kx_coeffs_roots;
    
end

%--
% might not be this simple -- kx root values are likely not 'ordered' (MS)
k1 = kx_arr(:,1);
k2 = kx_arr(:,2);
k3 = kx_arr(:,3);
k4 = kx_arr(:,4);

%%
%--
% plot kx's
figure(5)
plot(xax, real(k1),'.k')

hold on

plot(xax, imag(k1),'.r')
plot(xax, real(k3),'dk','MarkerSize',3)
plot(xax, imag(k3),'dr','MarkerSize',3)
legend('Re[k1]', 'Im[k1]', 'Re[k3]', 'Im[k3]')
% ylim([-2.5, 2.5])
xlabel('Position, x ($m$)')
% ylabel('Real(kx) (/m)')

hold off

figure(6)
plot(xax, real(k2),'.k')

hold on

plot(xax, imag(k2),'.r')
plot(xax, real(k4),'dk','MarkerSize',3)
plot(xax, imag(k4),'dr','MarkerSize',3)
legend('Re[k2]', 'Im[k2]', 'Re[k4]', 'Im[k4]')
% ylim([-10.0, 10.0])
xlabel('Position, x ($m$)')
% ylabel('Imag(kx) (/m)')

hold off






