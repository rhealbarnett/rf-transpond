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
    % the determinant of the above equation will be a quartic in kx -- the
    % dispersion relation
    kx_quart = det(wave_eq);
    
    %--
    % coeffs + 'All' finds the polynomial coeffients on the highest to
    % lowest order terms (ie for ax^4 + bx^3 ... etc they are ordered [a,
    % b, c, d, e]
    kx_coeffs = coeffs(kx_quart, 'All');
    
    %--
    % the roots function uses the polynomial coefficients, in order highest
    % to lowest, to determine the polynomial roots. 
    kx_coeffs_roots = roots(kx_coeffs);
    
    %--
    % store the four roots
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

imme = find(imag(me)==0);
imme = imme(end)+1;

figure(5)
semilogy(xax, real(k1),'.k')

hold on

semilogy(xax, imag(k1),'.r')
semilogy(xax, real(k3),'dk','MarkerSize',3)
semilogy(xax, imag(k3),'dr','MarkerSize',3)
legend('Re[k1]', 'Im[k1]', 'Re[k3]', 'Im[k3]')
% ylim([-2.5, 2.5])
xlabel('Position ($m$)','Fontsize',16)
% ylabel('Real(kx) (/m)')
vline(xax(imme),'--k')

hold off

figure(6)
semilogy(xax, real(k2),'.k')

hold on

semilogy(xax, imag(k2),'.r')
semilogy(xax, real(k4),'dk','MarkerSize',3)
semilogy(xax, imag(k4),'dr','MarkerSize',3)
legend('Re[k2]', 'Im[k2]', 'Re[k4]', 'Im[k4]')
% ylim([-10.0, 10.0])
xlabel('Position ($m$)','Fontsize',16)
% ylabel('Imag(kx) (/m)')
vline(xax(imme),'--k')

hold off

%%
%--
% transform data for log plot
y1 = sign(k1).*log10(abs(k1));
y2 = sign(k2).*log10(abs(k2));
y3 = sign(k3).*log10(abs(k3));
y4 = sign(k4).*log10(abs(k4));

figure(7)
plot(xax,real(y1),'.k')

hold on

plot(xax,imag(y1),'.r')
plot(xax,real(y3),'dk','MarkerSize',3)
plot(xax,imag(y3),'dr','MarkerSize',3)
legend('Re[k1]', 'Im[k1]', 'Re[k3]', 'Im[k3]')
xlabel('Position ($m$)','Fontsize',16)
vline(xax(imme),'--k')
% yticklabels({'-10$^{3}$','-10$^{2}$','-10$^{1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'})
ylabel('log$_{10}|$k$_x|$','Fontsize',16)

hold off

figure(8)
plot(xax,real(y2),'.k')

hold on

plot(xax,imag(y2),'.r')
plot(xax,real(y4),'dk','MarkerSize',3)
plot(xax,imag(y4),'dr','MarkerSize',3)
legend('Re[k2]', 'Im[k2]', 'Re[k4]', 'Im[k4]')
xlabel('Position ($m$)','Fontsize',16)
vline(xax(imme),'--k')
% yticklabels({'-10$^{3}$','-10$^{2}$','-10$^{1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'})
ylabel('log$_{10}|$k$_x|$','Fontsize',16)

hold off

