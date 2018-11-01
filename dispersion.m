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

% kpara = kz;
npara = c0*k_para/om;

r_arr = s_arr + d_arr;
l_arr = s_arr - d_arr;

a1 = s_arr;
b1 = r_arr.*l_arr + p_arr.*s_arr - npara^2*(p_arr + s_arr);
c1 = p_arr.*((npara^2 - r_arr).*(npara^2 - l_arr));

ns_p1 = (b1 - sqrt(b1.^2 - 4.0*a1.*c1))./(2.0*a1);
ns_p2 = (b1 + sqrt(b1.^2 - 4.0*a1.*c1))./(2.0*a1);

np11 = sqrt(ns_p1);
np12 = -sqrt(ns_p1);
np21 = sqrt(ns_p2);
np22 = -sqrt(ns_p2);

kp11 = np11*om/c0;
kp12 = np12*om/c0;
kp21 = np21*om/c0;
kp22 = np22*om/c0;

kperp_arr(:,1) = kp11;
kperp_arr(:,2) = kp12;
kperp_arr(:,3) = kp21;
kperp_arr(:,4) = kp22;

ns_s = -(npara.^2 - s_arr).*(p_arr./s_arr);
ns_f = -((npara.^2 - r_arr).*(npara.^2 - l_arr))./(npara.^2 - s_arr);
n_s1 = sqrt(ns_s); n_s2 = -sqrt(ns_s); 
n_f1 = sqrt(ns_f); n_f2 = -sqrt(ns_f);

k_s1 = n_s1*om/c0; k_s2 = n_s2*om/c0; 
k_f1 = n_f1*om/c0; k_f2 = n_f2*om/c0;

%%
%--
% plot kperps's

%--
% transform data for log plot
yp11 = sign(kp11).*log10(abs(kp11));
yp12 = sign(kp12).*log10(abs(kp12));
yp21 = sign(kp21).*log10(abs(kp21));
yp22 = sign(kp22).*log10(abs(kp22));

ks1 = sign(k_s1).*log10(abs(k_s1));
ks2 = sign(k_s2).*log10(abs(k_s2));
kf1 = sign(k_f1).*log10(abs(k_f1));
kf2 = sign(k_f2).*log10(abs(k_f2));

%%

figure(9)
plot(log10(n_new),real(yp11),'.k')

hold on

plot(log10(n_new),imag(yp11),'.r')
plot(log10(n_new),real(yp12),'dk','MarkerSize',3)
plot(log10(n_new),imag(yp12),'dr','MarkerSize',3)
legend('Re[k_{\perp1}]', 'Im[k_{\perp1}]', 'Re[k_{\perp1}]', 'Im[k_{\perp1}]')
xlabel('log_{10}|n|','Fontsize',16)
% vline(log10(N0(imme)),'--k')
ylabel('log_{10}|k_{\perp1}|','Fontsize',16)
xlim([log10(min(n_new)),log10(max(n_new))])
set(gca,'XDir','reverse');

hold off

figure(10)
plot(log10(n_new),real(yp21),'.k')

hold on

plot(log10(n_new),imag(yp21),'.r')
plot(log10(n_new),real(yp22),'dk','MarkerSize',3)
plot(log10(n_new),imag(yp22),'dr','MarkerSize',3)
legend('Re[k_{\perp2}]', 'Im[k_{\perp2}]', 'Re[k_{\perp2}]', 'Im[k_{\perp2}]')
xlabel('log_{10}|n|','Fontsize',16)
% vline(log10(N0(imme)),'--k') 
ylabel('log_{10}|k_{\perp2}|','Fontsize',16)
xlim([log10(min(n_new)),log10(max(n_new))])
set(gca,'XDir','reverse');

hold off

%%
% 
% figure(11)
% plot(log10(n_new),real(ks1),'.k')
% 
% hold on
% 
% plot(log10(n_new),imag(ks1),'.r')
% plot(log10(n_new),real(ks2),'dk','MarkerSize',3)
% plot(log10(n_new),imag(ks2),'dr','MarkerSize',3)
% legend('Re[k_{s1}]', 'Im[k_{s1}]', 'Re[k_{s1}]', 'Im[k_{s1}]')
% xlabel('log_{10}|n|','Fontsize',16)
% % vline(log10(N0(imme)),'--k')
% ylabel('log_{10}|k_{s1}|','Fontsize',16)
% xlim([log10(min(n_new)),log10(max(n_new))]);
% 
% hold off
% 
% figure(12)
% plot(log10(n_new),real(kf1),'.k')
% 
% hold on
% 
% plot(log10(n_new),imag(kf1),'.r')
% plot(log10(n_new),real(kf2),'dk','MarkerSize',3)
% plot(log10(n_new),imag(kf2),'dr','MarkerSize',3)
% legend('Re[k_{f2}]', 'Im[k_{f2}]', 'Re[k_{f2}]', 'Im[k_{f2}]')
% xlabel('log_{10}|n|','Fontsize',16)
% % vline(log10(N0(imme)),'--k') 
% ylabel('log_{10}|k_{f2}|','Fontsize',16)
% xlim([log10(min(n_new)),log10(max(n_new))]);
% 
% hold off
% 
% %%
% 
% figure(9)
% plot(log10(n_new),real(kp11),'.k')
% 
% hold on
% 
% plot(log10(n_new),imag(kp11),'.r')
% plot(log10(n_new),real(kp12),'dk','MarkerSize',3)
% plot(log10(n_new),imag(kp12),'dr','MarkerSize',3)
% legend('Re[k_{\perp1}]', 'Im[k_{\perp1}]', 'Re[k_{\perp1}]', 'Im[k_{\perp1}]')
% xlabel('log_{10}|n|','Fontsize',16)
% % vline(log10(N0(imme)),'--k')
% ylabel('log_{10}|k_{\perp1}|','Fontsize',16)
% % xlim([log10(min(n_new)),log10(max(n_new))]);
% 
% hold off
% 
% figure(10)
% plot(log10(n_new),real(kp21),'.k')
% 
% hold on
% 
% plot(log10(n_new),imag(kp21),'.r')
% plot(log10(n_new),real(kp22),'dk','MarkerSize',3)
% plot(log10(n_new),imag(kp22),'dr','MarkerSize',3)
% legend('Re[k_{\perp2}]', 'Im[k_{\perp2}]', 'Re[k_{\perp2}]', 'Im[k_{\perp2}]')
% xlabel('log_{10}|n|','Fontsize',16)
% % vline(log10(N0(imme)),'--k') 
% ylabel('log_{10}|k_{\perp2}|','Fontsize',16)
% % xlim([log10(min(n_new)),log10(max(n_new))]);
% 
% hold off
% 
% %%
% 
% % %%
% % 
% % figure(9)
% % plot(xax(npts-imme+1:imme-1),real(kp11(npts-imme+1:imme-1)),'.k')
% % 
% % hold on
% % 
% % plot(xax(npts-imme+1:imme-1),imag(kp11(npts-imme+1:imme-1)),'.r')
% % plot(xax(npts-imme+1:imme-1),real(kp12(npts-imme+1:imme-1)),'dk','MarkerSize',3)
% % plot(xax(npts-imme+1:imme-1),imag(kp12(npts-imme+1:imme-1)),'dr','MarkerSize',3)
% % legend('Re[k$_{\perp1}$]', 'Im[k$_{\perp1}$]', 'Re[k$_{\perp1}$]', 'Im[k$_{\perp1}$]')
% % xlabel('Position','Fontsize',16)
% % % vline(xax(imme),'--k')
% % ylabel('k$_{\perp1}$','Fontsize',16)
% % 
% % hold off
% % 
% % figure(10)
% % plot(xax(npts-imme+1:imme-1),real(kp21(npts-imme+1:imme-1)),'.k')
% % 
% % hold on
% % 
% % plot(xax(npts-imme+1:imme-1),imag(kp21(npts-imme+1:imme-1)),'.r')
% % plot(xax(npts-imme+1:imme-1),real(kp22(npts-imme+1:imme-1)),'dk','MarkerSize',3)
% % plot(xax(npts-imme+1:imme-1),imag(kp22(npts-imme+1:imme-1)),'dr','MarkerSize',3)
% % legend('Re[k$_{\perp2}$]', 'Im[k$_{\perp2}$]', 'Re[k$_{\perp2}$]', 'Im[k$_{\perp2}$]')
% % xlabel('Position','Fontsize',16)
% % % vline(xax(imme),'--k')
% % ylabel('k$_{\perp2}$','Fontsize',16)
% % 
% % hold off
% % 
% % 
% % 
% % 
