%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% values/parameters from van eester 2015  %
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

syms kx ky kz kperp kpara

%--
% initialise kx roots arrays, ensure they are complex
kperp_arr = zeros(npts, 4);
kperp_arr = complex(kperp_arr);
kp1_arr = zeros(npts, 1);
kp2_arr = zeros(npts, 1);

%-- 
% k matrix
a11 = -(ky^2 + kz^2);
a12 = -ky*kx;
a13 = kz*kx;
a21 = -ky*kx;
a22 = kz^2 + kx^2;
a23 = -ky*kz;
a31 = kz*kx;
a32 = -ky*kz;
a33 = -(ky^2 + kx^2);

a = [[a11, a12, a13]
    [a21, a22, a23]
    [a31, a32, a33]];

a = subs(a, [kx kz], [sqrt(kperp^2 - ky^2) kpara]);

kpara = 5.0;
ky = 0.0;

a = subs(a);

kperp_quart_arr = sym('K',[npts,1]);
check = zeros(npts, 4);
kperp_coeffs_arr = zeros(npts,5);

%-- 
% loop through density values

%--
% wave equation rhs
we_rhs = k0^2*cpdt;

for ii = 1:npts

    %--
    % set wave equation to zero to find determinant
    wave_eq = a - we_rhs(:,:,ii);

    %--
    % the determinant of the above equation will be a quartic in kx -- the
    % dispersion relation
    kperp_quart = det(wave_eq);
    kperp_quart_arr(ii,1) = kperp_quart;
    
    %--
    % coeffs + 'All' finds the polynomial coeffients on the highest to
    % lowest order terms (ie for ax^4 + bx^3 ... etc they are ordered [a,
    % b, c, d, e]
    kperp_coeffs = coeffs(kperp_quart, kperp, 'All');
    kperp_coeffs_arr(ii,:) = kperp_coeffs;
    
    c4 = kperp_coeffs(1);
    c2 = kperp_coeffs(3);
    c = kperp_coeffs(5);
    kperp1 = (-1.0*c2 - sqrt(c2^2 - 4.0*c4*c))/(2.0*c4);
    kperp2 = (-1.0*c2 + sqrt(c2^2 - 4.0*c4*c))/(2.0*c4);
    kp1_arr(ii,1) = kperp1;
    kp2_arr(ii,1) = kperp2;
    
    
    %--
    % the roots function uses the polynomial coefficients, in order highest
    % to lowest, to determine the polynomial roots. 
    kperp_coeffs_roots = roots(kperp_coeffs);
    
    %--
    % store the four roots
    kperp_arr(ii,:) = kperp_coeffs_roots;
    
    for kk=1:4
        check(ii,kk) = vpa(subs(kperp_quart,kperp,kperp_arr(ii,kk)));
    end
    
end

%--
% might not be this simple -- kx root values are likely not 'ordered' (MS)
k1 = kperp_arr(:,1);
k2 = kperp_arr(:,2);
k3 = kperp_arr(:,3);
k4 = kperp_arr(:,4);

%%
%--
% plot kx's
imme = find(imag(me)==0);
imme = imme(end)+1;

%--
% transform data for log plot
y1 = sign(k1).*log10(abs(k1));
y2 = sign(k2).*log10(abs(k2));
y3 = sign(k3).*log10(abs(k3));
y4 = sign(k4).*log10(abs(k4));

sw1 = sqrt(kp1_arr);
sw2 = -sqrt(kp1_arr);
fw1 = sqrt(kp2_arr);
fw2 = -sqrt(kp2_arr);
sf_arr = zeros(npts,4);
sf_arr(:,1) = sw1;
sf_arr(:,2) = sw2;
sf_arr(:,3) = fw1;
sf_arr(:,4) = fw2;


s1 = sign(sw1).*log10(abs(sw1));
s2 = sign(sw2).*log10(abs(sw1));
f1 = sign(fw1).*log10(abs(fw1));
f2 = sign(fw2).*log10(abs(fw1));


%%

figure(9)
plot(log10(N0),real(s1),'.k')

hold on

plot(log10(N0),imag(s1),'.r')
plot(log10(N0),real(s2),'dk','MarkerSize',3)
plot(log10(N0),imag(s2),'dr','MarkerSize',3)
legend('Re[k$_{\perp1}$]', 'Im[k$_{\perp1}$]', 'Re[k$_{\perp1}$]', 'Im[k$_{\perp1}$]')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
vline(log10(N0(imme)),'--k')
ylabel('log$_{10}|$k$_{\perp1}|$','Fontsize',16)

hold off

figure(15)
plot(log10(N0),real(f1),'.k')

hold on

plot(log10(N0),imag(f1),'.r')
plot(log10(N0),real(f2),'dk','MarkerSize',3)
plot(log10(N0),imag(f2),'dr','MarkerSize',3)
legend('Re[k$_{\perp2}$]', 'Im[k$_{\perp2}$]', 'Re[k$_{\perp2}$]', 'Im[k$_{\perp2}$]')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
vline(log10(N0(imme)),'--k')
ylabel('log$_{10}|$k$_{\perp2}|$','Fontsize',16)

hold off

