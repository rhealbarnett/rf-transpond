%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% values/parameters from van eester 2015  %
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

perp = 1;
para = 0;

r_arr = s_arr + d_arr;
l_arr = s_arr - d_arr;

kp1_arr = zeros(npts, 1);
kp2_arr = zeros(npts, 1);

k_para11 = zeros(npts,length(n_para));
k_para12 = zeros(npts,length(n_para));
k_para21 = zeros(npts,length(n_para));
k_para22 = zeros(npts,length(n_para));

if perp
    
    syms kperp

    %--
    % initialise kx roots arrays, ensure they are complex
    kperp_arr = zeros(npts, 4);
    kperp_arr = complex(kperp_arr);

    for ii=1:length(n_para)
        a1 = s_arr;
        b1 = r_arr.*l_arr + p_arr.*s_arr - n_para(1,ii)^2*(p_arr + s_arr);
        c1 = p_arr.*((n_para(1,ii)^2 - r_arr).*(n_para(1,ii)^2 - l_arr));

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
        
        kpara11(:,ii) = kp11;
        kpara12(:,ii) = kp12;
        kpara21(:,ii) = kp21;
        kpara22(:,ii) = kp22;
% 
%         ns_s = -(npara.^2 - s_arr).*(p_arr./s_arr);
%         ns_f = -((npara.^2 - r_arr).*(npara.^2 - l_arr))./(npara.^2 - s_arr);
%         n_s1 = sqrt(ns_s); n_s2 = -sqrt(ns_s); 
%         n_f1 = sqrt(ns_f); n_f2 = -sqrt(ns_f);
% 
%         k_s1 = n_s1*om/c0; k_s2 = n_s2*om/c0; 
%         k_f1 = n_f1*om/c0; k_f2 = n_f2*om/c0;
% 
%         ks1 = sign(k_s1).*log10(abs(k_s1));
%         ks2 = sign(k_s2).*log10(abs(k_s2));
%         kf1 = sign(k_f1).*log10(abs(k_f1));
%         kf2 = sign(k_f2).*log10(abs(k_f2));
    end
    
    nsize = size(n_new);
    
    if nsize(1)~=1
        n_new = n_new(1,:);
    else
    end
        
    if length(n_para)==1
        figure(1)
        subplot(1,2,1)
        semilogx(n_new, real(kp21),'.k')
        hold on
        semilogx(n_new, real(kp22),'.k')
        semilogx(n_new, imag(kp21),'.r')
        semilogx(n_new, imag(kp22),'.r')
        xlabel('n (m^{-3})')
        ylabel('k_{\perp} (m^{-1})')

        subplot(1,2,2)
        semilogx(n_new, real(kp11),'.k')
        hold on
        semilogx(n_new, real(kp12),'.k')
        semilogx(n_new, imag(kp11),'.r')
        semilogx(n_new, imag(kp12),'.r')
        xlabel('n (m^{-3})')
        ylabel('k_{\perp} (m^{-1})')
        
    elseif length(n_para)~=1
        
        figure(2)
        subplot(2,2,1)
        contourf(log10(n_new),(k_para),real(kpara11)','Linecolor','none')
        rc11=colorbar;
        set(gca,'xtick',[])
        title('Re[k_{\perp 11}]')
        ylabel('k_{||} (m^{-1})')

        subplot(2,2,2)
        contourf(log10(n_new),(k_para),real(kpara12)','Linecolor','none')
        rc11 = colorbar;
        title('Re[k_{\perp 12}]')
        ylabel('k_{||} (m^{-1})')
        xlabel('log_{10}n')

        subplot(2,2,3)
        contourf(log10(n_new),(k_para),real(kpara21)','Linecolor','none')
        rc11 = colorbar;
        title('Re[k_{\perp 21}]')
        ylabel('k_{||} (m^{-1})')
        xlabel('log_{10}n')

        subplot(2,2,4)
        contourf(log10(n_new),(k_para),real(kpara22)','Linecolor','none')
        rc11 = colorbar;
        title('Re[k_{\perp 22}]')
        ylabel('k_{||} (m^{-1})')
        xlabel('log_{10}n')
        
    end
    
elseif para
    
    syms kpara

    kpara_arr = zeros(npts, 4, npts);
    kpara_arr = complex(kpara_arr);
    
    %--
    % initialise kx roots arrays, ensure they are complex
    for ii=1:length(n_perp)

        a1 = p_arr;
        b1 = (2.0*p_arr.*s_arr - n_perp(1,ii)^2*(p_arr + s_arr));
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
    
    if length(n_perp)==1
        figure(2)
        subplot(1,2,1)
        semilogx(n_new, real(kp21),'.k')
        hold on
        semilogx(n_new, real(kp22),'.k')
        semilogx(n_new, imag(kp21),'.r')
        semilogx(n_new, imag(kp22),'.r')
        xlabel('n (m^{-3})')
        ylabel('k_{||} (m^{-1})')

        subplot(1,2,2)
        semilogx(n_new, real(kp11),'.k')
        hold on
        semilogx(n_new, real(kp12),'.k')
        semilogx(n_new, imag(kp11),'.r')
        semilogx(n_new, imag(kp12),'.r')
        xlabel('n (m^{-3})')
        ylabel('k_{||} (m^{-1})')
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
plot(log10(Nmax)*ones(1,npts),log10(linspace(0.01,100,npts)),'b','Linewidth',3)
plot(log10(1.0e17)*ones(1,npts),log10(linspace(0.01,100,npts)),'b','Linewidth',3)
if perp
    legend('Re[k_{\perp1}]', 'Im[k_{\perp1}]', 'Re[k_{\perp1}]', 'Im[k_{\perp1}]')
    ylabel('log_{10}|k_{\perp1}|','Fontsize',16)
elseif para
    legend('Re[k_{||1}]', 'Im[k_{||1}]', 'Re[k_{||1}]', 'Im[k_{||1}]')
    ylabel('log_{10}|k_{||1}|','Fontsize',16)
end
xlabel('log_{10}|n|','Fontsize',16)
% vline(log10(Nmin),'--b')
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

figure(3)
subplot(2,2,1)
contourf(log10(n_new),(k_perp),real(kpara11)',levels,'Linecolor','none')
colormap(gca,flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(rc11,'Re[k_{||11}] (m^{-1})')
set(gca,'xtick',[])
title('Re[k_{||11}]')
ylabel('k_{\perp} (m^{-1})')
rc11.Ticks = linspace(0,14,5);


subplot(2,2,2)
contourf(log10(n_new),(k_perp),real(kpara12)',levels,'Linecolor','none')
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
contourf(log10(n_new),(k_perp),real(kpara21)',levels,'Linecolor','none')
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
contourf(log10(n_new),(k_perp),real(kpara22)',levels,'Linecolor','none')
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

levims11 = linspace(-max(abs(imag(kpara11(:)))),max(abs(imag(kpara11(:)))),50);
levims22 = linspace(-max(abs(imag(kpara22(:)))),max(abs(imag(kpara22(:)))),50);

figure(4)
subplot(2,2,1)
contourf(log10(n_new),(k_perp),imag(kpara11)',levims11,'Linecolor','none')
colormap(gca,flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(rc11,'Re[k_{||11}] (m^{-1})')
set(gca,'xtick',[])
title('Im[k_{||11}]')
ylabel('k_{\perp} (m^{-1})')
% rc11.Ticks = linspace(0,14,5);


subplot(2,2,2)
contourf(log10(n_new),(k_perp),imag(kpara12)',levims11,'Linecolor','none')
colormap(gca,flipud(c(41:64,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Im[k_{||12}]')
ylabel('k_{\perp} (m^{-1})')
% xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

subplot(2,2,3)
contourf(log10(n_new),(k_perp),imag(kpara21)',levims22,'Linecolor','none')
colormap(gca,flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Im[k_{||21}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

subplot(2,2,4)
contourf(log10(n_new),(k_perp),imag(kpara22)',levims22,'Linecolor','none')
colormap(gca,flipud(c(41:64,:)));
rc11 = colorbar;
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Im[k_{||22}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

%%

c = jet(64);

levels11 = linspace(-max(real(kpara11(:))),max(real(kpara11(:))),50);
levims11 = linspace(-max(abs(imag(kpara11(:)))),max(abs(imag(kpara11(:)))),50);
levels22 = linspace(-max(real(kpara21(:))),max(real(kpara21(:))),500);
levims22 = linspace(-max(abs(imag(kpara21(:)))),max(abs(imag(kpara21(:)))),50);

figure(5)
subplot(2,2,1)
contourf(log10(n_new),(k_perp),real(kpara11)',levels11,'Linecolor','none')
colormap(gca,flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(rc11,'Re[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
title('Re[k_{||11}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')

subplot(2,2,3)
contourf(log10(n_new),(k_perp),imag(kpara11)',levims11,'Linecolor','none')
colormap(gca,flipud(c(1:40,:)));
rc11 = colorbar;
% ylabel(rc11,'Re[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
title('Im[k_{||11}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,14,5);

subplot(2,2,2)
contourf(log10(n_new),(k_perp),real(kpara21)',levels22,'Linecolor','none')
if real(kpara21(:))==0
    shading flat
    colormap(gca,flipud(c(1:40,:)));
    rc11 = colorbar;
    lim = caxis;
    caxis([0 1]);
else
    colormap(gca,flipud(c(1:40,:)));
    rc11 = colorbar;
end
%flipud(c(1:40,:)));
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Re[k_{||21}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

subplot(2,2,4)
contourf(log10(n_new),(k_perp),imag(kpara21)',levims22,'Linecolor','none')
colormap(gca,flipud(c(1:40,:)));
if imag(kpara21(:))==0
    shading flat
    colormap(gca,flipud(c(1:40,:)));
    rc11 = colorbar;
    lim = caxis;
    caxis([0 1]);
else
    colormap(gca,flipud(c(1:40,:)));
    rc11 = colorbar;
end
% ylabel(ic11,'Im[k_{||11}] (m^{-1})')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Im[k_{||21}]')
ylabel('k_{\perp} (m^{-1})')
xlabel('log_{10}n')
% rc11.Ticks = linspace(0,3,4);

