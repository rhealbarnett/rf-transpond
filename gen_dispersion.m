%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% can solve when B0 is not aligned with z %
% rlbarnett c3149416, 010218              %
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

% kx_quart_arr = sym('K',[npts,1]);
check = zeros(npts, 4);

%--
% wave equation rhs
we_rhs = k0^2*cpdt;

%-- 
% loop through density values

for ii = 1:npts

    %--
    % set wave equation to zero to find determinant
    wave_eq = a - we_rhs(:,:,ii);

    %--
    % the determinant of the above equation will be a quartic in kx -- the
    % dispersion relation
    kx_quart = det(wave_eq);
%     kx_quart_arr(ii,1) = kx_quart;
    
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
    
    for kk=1:4
        check(ii,kk) = vpa(subs(kx_quart,kx,kx_arr(ii,kk)));
    end
    
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


%%
%--
% transform data for log plot
y1 = sign(k1).*log10(abs(k1));
y2 = sign(k2).*log10(abs(k2));
y3 = sign(k3).*log10(abs(k3));
y4 = sign(k4).*log10(abs(k4));

figure(7)
plot(log10(N0),real(y1),'.k')

hold on

plot(log10(N0),imag(y1),'.r')
plot(log10(N0),real(y3),'dk','MarkerSize',3)
plot(log10(N0),imag(y3),'dr','MarkerSize',3)
legend('Re[k1]', 'Im[k1]', 'Re[k3]', 'Im[k3]')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
vline(log10(N0(imme)),'--k')
% yticklabels({'-10$^{3}$','-10$^{2}$','-10$^{1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'})
ylabel('log$_{10}|$k$_x|$','Fontsize',16)

hold off

figure(8)
plot(log10(N0),real(y2),'.k')

hold on

plot(log10(N0),imag(y2),'.r')
plot(log10(N0),real(y4),'dk','MarkerSize',3)
plot(log10(N0),imag(y4),'dr','MarkerSize',3)
legend('Re[k2]', 'Im[k2]', 'Re[k4]', 'Im[k4]')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
vline(log10(N0(imme)),'--k')
% yticklabels({'-10$^{3}$','-10$^{2}$','-10$^{1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'})
ylabel('log$_{10}|$k$_x|$','Fontsize',16)

hold off