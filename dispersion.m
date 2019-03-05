%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% values/parameters from van eester 2015  %
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

perp = 0;
para = 1;

r_arr = s_arr + d_arr;
l_arr = s_arr - d_arr;

kp1_arr = zeros(npts, 1);
kp2_arr = zeros(npts, 1);


if perp
    
    syms kperp

    %--
    % initialise kx roots arrays, ensure they are complex
    kperp_arr = zeros(npts, 4);
    kperp_arr = complex(kperp_arr);

    a1 = s_arr;
    b1 = r_arr.*l_arr + p_arr.*s_arr - npara^2*(p_arr + s_arr);
    c1 = p_arr.*((npara^2 - r_arr).*(npara^2 - l_arr));
    
    ns_s = -(npara.^2 - s_arr).*(p_arr./s_arr);
    ns_f = -((npara.^2 - r_arr).*(npara.^2 - l_arr))./(npara.^2 - s_arr);
    n_s1 = sqrt(ns_s); n_s2 = -sqrt(ns_s); 
    n_f1 = sqrt(ns_f); n_f2 = -sqrt(ns_f);

    k_s1 = n_s1*om/c0; k_s2 = n_s2*om/c0; 
    k_f1 = n_f1*om/c0; k_f2 = n_f2*om/c0;
    
    ks1 = sign(k_s1).*log10(abs(k_s1));
    ks2 = sign(k_s2).*log10(abs(k_s2));
    kf1 = sign(k_f1).*log10(abs(k_f1));
    kf2 = sign(k_f2).*log10(abs(k_f2));
    
elseif para
    
    syms kpara

    kpara_arr = zeros(npts, 4, npts);
    kpara_arr = complex(kpara_arr);
    
    %--
    % initialise kx roots arrays, ensure they are complex
    for ii=1:npts

        a1 = p_arr;
        b1 = 2.0*p_arr.*s_arr - n_perp(1,ii)^2*(p_arr + s_arr);
        c1 = (n_perp(1,ii)^2 - p_arr).*(s_arr*n_perp(1,ii)^2 - r_arr.*l_arr);
        
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
        
        kpara_arr(:,1,ii) = kp11;
        kpara_arr(:,2,ii) = kp12;
        kpara_arr(:,3,ii) = kp21;
        kpara_arr(:,4,ii) = kp22;
        
        kpara11(:,ii) = kp11;
        kpara12(:,ii) = kp12;
        kpara21(:,ii) = kp21;
        kpara22(:,ii) = kp22;
        
    end
    
end


kperp_arr(:,1) = kp11;
kperp_arr(:,2) = kp12;
kperp_arr(:,3) = kp21;
kperp_arr(:,4) = kp22;


%%
%--
% plot k's

%--
% transform data for log plot
yp11 = sign(kp11).*log10(abs(kp11));
yp12 = sign(kp12).*log10(abs(kp12));
yp21 = sign(kp21).*log10(abs(kp21));
yp22 = sign(kp22).*log10(abs(kp22));


%%

figure(9)
plot(log10(n_new),real(yp11),'.k','Markersize',20)

hold on

plot(log10(n_new),imag(yp11),'.r','Markersize',20)
plot(log10(n_new),real(yp12),'dk','MarkerSize',7)
plot(log10(n_new),imag(yp12),'dr','MarkerSize',7)
if perp
    legend('Re[k_{\perp1}]', 'Im[k_{\perp1}]', 'Re[k_{\perp1}]', 'Im[k_{\perp1}]')
    ylabel('log_{10}|k_{\perp1}|','Fontsize',16)
elseif para
    legend('Re[k_{||1}]', 'Im[k_{||1}]', 'Re[k_{||1}]', 'Im[k_{||1}]')
    ylabel('log_{10}|k_{||1}|','Fontsize',16)
end
xlabel('log_{10}|n|','Fontsize',16)
% vline(log10(N0(imme)),'--k')
xlim([log10(min(n_new)),log10(max(n_new))])
set(gca,'XDir','reverse');

hold off

figure(10)
plot(log10(n_new),real(yp21),'.k','Markersize',20)

hold on

plot(log10(n_new),imag(yp21),'.r','Markersize',20)
plot(log10(n_new),real(yp22),'dk','MarkerSize',7)
plot(log10(n_new),imag(yp22),'dr','MarkerSize',7)
if perp
    legend('Re[k_{\perp2}]', 'Im[k_{\perp2}]', 'Re[k_{\perp2}]', 'Im[k_{\perp2}]')
    ylabel('log_{10}|k_{\perp2}|','Fontsize',16)
elseif para
    legend('Re[k_{||2}]', 'Im[k_{||2}]', 'Re[k_{||2}]', 'Im[k_{||2}]')
    ylabel('log_{10}|k_{||2}|','Fontsize',16)
end
xlabel('log_{10}|n|','Fontsize',16)
% vline(log10(N0(imme)),'--k') 
xlim([log10(min(n_new)),log10(max(n_new))])
set(gca,'XDir','reverse');

hold off

%%

c = jet(64);

levels = linspace(-max(real(kpara11(:))),max(real(kpara11(:))),50);
levims = linspace(-max(imag(kpara11(:))),max(imag(kpara11(:))),50);

subplot(2,2,1)
contourf(log10(n_new),(k_perp),real(kpara11),levels,'Linecolor','none')
colormap(gca,flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(rc11,'Re[k_{||11}] (m^{-1})')
set(gca,'xtick',[])
title('Re[k_{||11}]')
ylabel('k_{\perp} (m^{-1})')
rc11.Ticks = linspace(0,14,5);


subplot(2,2,2)
contourf(log10(n_new),(k_perp),real(kpara12),levels,'Linecolor','none')
colormap(gca,flipud(c(41:64,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Re[k_{||12}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

subplot(2,2,3)
contourf(log10(n_new),(k_perp),real(kpara21),levels,'Linecolor','none')
colormap(gca,c);%flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Re[k_{||21}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

subplot(2,2,4)
contourf(log10(n_new),(k_perp),real(kpara22),levels,'Linecolor','none')
colormap(gca,c);%flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Re[k_{||22}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

%%

c = jet(64);

levels = linspace(-max(real(kpara11(:))),max(real(kpara11(:))),50);
levims = linspace(-max(imag(kpara11(:))),max(imag(kpara11(:))),50);

figure(2)
subplot(2,2,1)
contourf(log10(n_new),(k_perp),imag(kpara11),levels,'Linecolor','none')
colormap(gca,flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(rc11,'Re[k_{||11}] (m^{-1})')
set(gca,'xtick',[])
title('Im[k_{||11}]')
ylabel('k_{\perp} (m^{-1})')
rc11.Ticks = linspace(0,14,5);


subplot(2,2,2)
contourf(log10(n_new),(k_perp),imag(kpara12),levels,'Linecolor','none')
colormap(gca,flipud(c(41:64,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Im[k_{||12}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

subplot(2,2,3)
contourf(log10(n_new),(k_perp),imag(kpara21),levels,'Linecolor','none')
colormap(gca,c);%flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Im[k_{||21}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

subplot(2,2,4)
contourf(log10(n_new),(k_perp),imag(kpara22),levels,'Linecolor','none')
colormap(gca,c);%flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Im[k_{||22}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

%%

% subplot(2,2,3)
% contourf(log10(n_new),(k_perp),real(kpara12),levels,'Linecolor','none')
% colormap(gca,flipud(c(41:64,:)));
% rc12 = colorbar;
% ylabel(rc12,'Re[k_{||12}] (m^{-1})')
% xlabel('log_{10} n')
% ylabel('k_{\perp} (m^{-1})')
% 
% subplot(2,2,4)
% contourf(log10(n_new),(k_perp),imag(kpara12),levims,'Linecolor','none')
% colormap(gca,flipud(c(41:64,:)));
% ic12 = colorbar;
% ylabel(ic12,'Im[k_{||12}] (m^{-1})')
% % set(gca,'ytick',[])
% xlabel('log_{10} n')



% %%
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

%%
% 
% figure(11)
% plot(nxax,real(kp11),'.k')
% 
% hold on
% 
% plot(nxax,imag(kp11),'.r')
% plot(nxax,real(kp12),'dk','MarkerSize',3)
% plot(nxax,imag(kp12),'dr','MarkerSize',3)
% legend('Re[k_{\perp1}]', 'Im[k_{\perp1}]', 'Re[k_{\perp1}]', 'Im[k_{\perp1}]')
% xlabel('Position (m)','Fontsize',16)
% % vline(log10(N0(imme)),'--k')
% ylabel('k_{\perp1}','Fontsize',16)
% xlim([xmin,xmax]);
% % set(gca,'XDir','reverse');
% 
% hold off
% 
% figure(12)
% plot(nxax,real(kp21),'.k')
% 
% hold on
% 
% plot(nxax,imag(kp21),'.r')
% plot(nxax,real(kp22),'dk','MarkerSize',3)
% plot(nxax,imag(kp22),'dr','MarkerSize',3)
% legend('Re[k_{\perp2}]', 'Im[k_{\perp2}]', 'Re[k_{\perp2}]', 'Im[k_{\perp2}]')
% xlabel('Position (m)','Fontsize',16)
% % vline(log10(N0(imme)),'--k') 
% ylabel('k_{\perp2}','Fontsize',16)
% xlim([xmin,xmax]);
% % set(gca,'XDir','reverse');
% 
% hold off

%%
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
