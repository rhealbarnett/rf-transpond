%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% values/parameters from van eester 2015  %
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

syms kperp

%--
% initialise kx roots arrays, ensure they are complex
kperp_arr = zeros(npts, 4);
kperp_arr = complex(kperp_arr);
kp1_arr = zeros(npts, 1);
kp2_arr = zeros(npts, 1);

kpara = kz;

a1 = s_arr;
b1 = r_arr.*l_arr + p_arr.*s_arr - kpara^2*(p_arr + s_arr);
c1 = p_arr.*((kpara^2 - r_arr).*(kpara^2 - l_arr));

ks_p1 = (b1 - sqrt(b1.^2 - 4.0*a1.*c1))./(2.0*a1);
ks_p2 = (b1 + sqrt(b1.^2 - 4.0*a1.*c1))./(2.0*a1);

kp11 = sqrt(ks_p1);
kp12 = -sqrt(ks_p1);
kp21 = sqrt(ks_p2);
kp22 = -sqrt(ks_p2);

%%
%--
% plot kperps's
imme = find(imag(me)==0);
imme = imme(end)+1;

%--
% transform data for log plot
yp11 = sign(kp11).*log10(abs(kp11));
yp12 = sign(kp12).*log10(abs(kp12));
yp21 = sign(kp21).*log10(abs(kp21));
yp22 = sign(kp22).*log10(abs(kp22));

%%

figure(9)
plot(log10(N0),real(yp11),'.k')

hold on

plot(log10(N0),imag(yp11),'.r')
plot(log10(N0),real(yp12),'dk','MarkerSize',3)
plot(log10(N0),imag(yp12),'dr','MarkerSize',3)
legend('Re[k$_{\perp1}$]', 'Im[k$_{\perp1}$]', 'Re[k$_{\perp1}$]', 'Im[k$_{\perp1}$]')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
vline(log10(N0(imme)),'--k')
ylabel('log$_{10}|$k$_{\perp1}|$','Fontsize',16)

hold off

figure(10)
plot(log10(N0),real(yp21),'.k')

hold on

plot(log10(N0),imag(yp21),'.r')
plot(log10(N0),real(yp22),'dk','MarkerSize',3)
plot(log10(N0),imag(yp22),'dr','MarkerSize',3)
legend('Re[k$_{\perp2}$]', 'Im[k$_{\perp2}$]', 'Re[k$_{\perp2}$]', 'Im[k$_{\perp2}$]')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
vline(log10(N0(imme)),'--k')
ylabel('log$_{10}|$k$_{\perp2}|$','Fontsize',16)

hold off


