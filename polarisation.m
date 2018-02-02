%-----------------------------------------%
% Calculate analytical wave polarisation  %
% rlbarnett c3149416, 010218              %
%-----------------------------------------%

exey_fast = -1i*(d_arr./(kpara^2 - s_arr.^2));

ezey_fast11 = -1i*kpara*kp11.*(d_arr./(p_arr.*(kpara^2 - s_arr)));
ezey_fast12 = -1i*kpara*kp12.*(d_arr./(p_arr.*(kpara^2 - s_arr)));
ezey_fast21 = -1i*kpara*kp21.*(d_arr./(p_arr.*(kpara^2 - s_arr)));
ezey_fast22 = -1i*kpara*kp22.*(d_arr./(p_arr.*(kpara^2 - s_arr)));

exez_slow11 = (kp11.*kpara)./(kpara^2 - s_arr);
exez_slow12 = (kp12.*kpara)./(kpara^2 - s_arr);
exez_slow21 = (kp21.*kpara)./(kpara^2 - s_arr);
exez_slow22 = (kp22.*kpara)./(kpara^2 - s_arr);

eyez_slow11 = 1i*(kpara*d_arr)./(kp11.*(kpara^2 - s_arr));
eyez_slow12 = 1i*(kpara*d_arr)./(kp12.*(kpara^2 - s_arr));
eyez_slow21 = 1i*(kpara*d_arr)./(kp21.*(kpara^2 - s_arr));
eyez_slow22 = 1i*(kpara*d_arr)./(kp22.*(kpara^2 - s_arr));

%%
%------------------------------------------------%
% set Ex=1, calculate Ey and Ez from ratios      %
%------------------------------------------------%

ey_fast = (kpara^2 - s_arr)./(-1i*d_arr);

ez_fast11 = -1i*kpara*kp11.*(d_arr./(p_arr.*(kpara^2 - s_arr))).*ey_fast;
ez_fast12 = -1i*kpara*kp12.*(d_arr./(p_arr.*(kpara^2 - s_arr))).*ey_fast;
ez_fast21 = -1i*kpara*kp21.*(d_arr./(p_arr.*(kpara^2 - s_arr))).*ey_fast;
ez_fast22 = -1i*kpara*kp22.*(d_arr./(p_arr.*(kpara^2 - s_arr))).*ey_fast;

ez_slow11 = (kpara^2 - s_arr)./(kp11.*kpara);
ez_slow12 = (kpara^2 - s_arr)./(kp12.*kpara);
ez_slow21 = (kpara^2 - s_arr)./(kp21.*kpara);
ez_slow22 = (kpara^2 - s_arr)./(kp22.*kpara);

ey_slow11 = 1i*((kpara*d_arr)./(kp11.*(kpara^2 - s_arr))).*ez_slow11;
ey_slow12 = 1i*((kpara*d_arr)./(kp12.*(kpara^2 - s_arr))).*ez_slow12;
ey_slow21 = 1i*((kpara*d_arr)./(kp21.*(kpara^2 - s_arr))).*ez_slow21;
ey_slow22 = 1i*((kpara*d_arr)./(kp22.*(kpara^2 - s_arr))).*ez_slow22;

%%

figure(20)
set(gcf,'Position',[4 68 857 886])
suptitle('Fast wave polarisation')
subplot(3,1,1)
plot(log10(N0(1:imme)),real(exey_fast(1:imme)),'.k')
hold on

plot(log10(N0(1:imme)),imag(exey_fast(1:imme)),'.r')
xlim([min(log10(N0)) log10(N0(imme))])
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
ylabel('E$_x$/E$_y$','Fontsize',16)
legend('Re[E$_x$/E$_y$]','Im[E$_x$/E$_y$]')
hold off

subplot(3,1,2)
plot(log10(N0(1:imme)),real(ezey_fast11(1:imme)),'.k')
hold on

plot(log10(N0(1:imme)),imag(ezey_fast11(1:imme)),'.r')
plot(log10(N0(1:imme)),real(ezey_fast12(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(ezey_fast12(1:imme)),'dr','MarkerSize',3)
xlim([min(log10(N0)) log10(N0(imme))])
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
ylabel('E$_z$/E$_y$','Fontsize',16)
legend('Re[E$_z$/E$_y$]$_{\perp1,1}$','Im[E$_z$/E$_y$]$_{\perp1,1}$','Re[E$_z$/E$_y$]$_{\perp1,2}$','Im[E$_z$/E$_y$]$_{\perp1,2}$')
hold off

subplot(3,1,3)
plot(log10(N0(1:imme)),real(ezey_fast21(1:imme)),'.k')
hold on

plot(log10(N0(1:imme)),imag(ezey_fast21(1:imme)),'.r')
plot(log10(N0(1:imme)),real(ezey_fast22(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(ezey_fast22(1:imme)),'dr','MarkerSize',3)
xlim([min(log10(N0)) log10(N0(imme))])
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
ylabel('E$_z$/E$_y$','Fontsize',16)
legend('Re[E$_z$/E$_y$]$_{\perp2,1}$','Im[E$_z$/E$_y$]$_{\perp2,1}$','Re[E$_z$/E$_y$]$_{\perp2,2}$','Im[E$_z$/E$_y$]$_{\perp2,2}$')
hold off

figure(21)
set(gcf,'Position',[865 69 813 886])
suptitle('Slow wave polarisation')
subplot(4,1,1)
plot(log10(N0(1:imme)),real(exez_slow11(1:imme)),'.k')
hold on

plot(log10(N0(1:imme)),imag(exez_slow11(1:imme)),'.r')
plot(log10(N0(1:imme)),real(exez_slow12(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(exez_slow12(1:imme)),'dr','MarkerSize',3)
xlim([min(log10(N0)) log10(N0(imme))])
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
ylabel('E$_x$/E$_z$','Fontsize',16)
legend('Re[E$_x$/E$_z$]$_{\perp1,1}$','Im[E$_x$/E$_z$]$_{\perp1,1}$','Re[E$_x$/E$_z$]$_{\perp1,2}$','Im[E$_x$/E$_z$]$_{\perp1,2}$')
hold off

subplot(4,1,2)
plot(log10(N0(1:imme)),real(exez_slow21(1:imme)),'.k')
hold on

plot(log10(N0(1:imme)),imag(exez_slow21(1:imme)),'.r')
plot(log10(N0(1:imme)),real(exez_slow22(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(exez_slow22(1:imme)),'dr','MarkerSize',3)
xlim([min(log10(N0)) log10(N0(imme))])
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
ylabel('E$_x$/E$_z$','Fontsize',16)
legend('Re[E$_x$/E$_z$]$_{\perp2,1}$','Im[E$_x$/E$_z$]$_{\perp2,1}$','Re[E$_x$/E$_z$]$_{\perp2,2}$','Im[E$_x$/E$_z$]$_{\perp2,2}$')
hold off

subplot(4,1,3)
plot(log10(N0(1:imme)),real(eyez_slow11(1:imme)),'.k')
hold on

plot(log10(N0(1:imme)),imag(eyez_slow11(1:imme)),'.r')
plot(log10(N0(1:imme)),real(eyez_slow12(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(eyez_slow12(1:imme)),'dr','MarkerSize',3)
xlim([min(log10(N0)) log10(N0(imme))])
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
ylabel('E$_y$/E$_z$','Fontsize',16)
legend('Re[E$_y$/E$_z$]$_{\perp1,1}$','Im[E$_y$/E$_z$]$_{\perp1,1}$','Re[E$_y$/E$_z$]$_{\perp1,2}$','Im[E$_y$/E$_z$]$_{\perp1,2}$')
hold off

subplot(4,1,4)
plot(log10(N0(1:imme)),real(eyez_slow21(1:imme)),'.k')
hold on

plot(log10(N0(1:imme)),imag(eyez_slow21(1:imme)),'.r')
plot(log10(N0(1:imme)),real(eyez_slow22(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(eyez_slow22(1:imme)),'dr','MarkerSize',3)
xlim([min(log10(N0)) log10(N0(imme))])
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
ylabel('E$_y$/E$_z$','Fontsize',16)
legend('Re[E$_y$/E$_z$]$_{\perp2,1}$','Im[E$_y$/E$_z$]$_{\perp2,1}$','Re[E$_y$/E$_z$]$_{\perp2,2}$','Im[E$_y$/E$_z$]$_{\perp2,2}$')
hold off


