%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% values/parameters from van eester 2015  %
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

perp = 0;
para = 1;
plots = 1;

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
        
    if length(n_para)==1 && plots
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
        
    elseif length(n_para)~=1 && plots
         
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
        set(gca,'Fontsize',20)
        
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
    
    if length(n_perp)==1 && plots
        
        x0 = 0;
        y0 = 0;
        width = 1000;
        height = 500;
        
        figure(2)
        
        set(gcf,'Position',[x0 y0 width height],'Color','w')
        subplot(1,2,1)
        semilogx(n_new, real(kp21),'.k')
        hold on
        semilogx(n_new, real(kp22),'.k')
        semilogx(n_new, imag(kp21),'.r')
        semilogx(n_new, imag(kp22),'.r')
        semilogx((1.0e19)*ones(1,npts),(linspace(min(imag(kp22)),max(imag(kp21)),npts)),...
            'b--','Linewidth',1.5)
        semilogx((1.0e18)*ones(1,npts),(linspace(min(imag(kp22)),max(imag(kp21)),npts)),...
            'b--','Linewidth',1.5)
        ylim([-40 40])
        xlabel('$n$ (m$^{-3}$)','Interpreter','latex')
        ylabel('$k_{z}$ (m$^{-1}$)','Interpreter','latex')
        set(gca,'Fontsize',25)

        subplot(1,2,2)
        semilogx(n_new, real(kp11),'.k')
        hold on
        semilogx(n_new, real(kp12),'.k')
        semilogx(n_new, imag(kp11),'.r')
        semilogx(n_new, imag(kp12),'.r')
        semilogx((1.0e19)*ones(1,npts),(linspace(-60,60,npts)),...
            '--b','Linewidth',1.5)
        semilogx((1.0e18)*ones(1,npts),(linspace(-60,60,npts)),...
            '--b','Linewidth',1.5)
        ylim([-40 40])
        xlabel('$n$ (m$^{-3}$)','Interpreter','latex')
%         ylabel('$k_{z}$ (m$^{-1}$)','Interpreter','latex')
        yticks([])
        
        set(gca,'Fontsize',25)
        
%         export_fig('/Volumes/DATA/thesis/RFT/figs/kvsn_dispersion_lines.png',...
%             '-r300')
        
    elseif length(n_perp)~=1 && plots
        
        width = 1000;
        height = 800;
        x0 = 0;
        y0 = 0;
        
        indn_np = find(n_new<=1.0e17,1,'first');
        indky_np = find(ky<=30,1,'last');
        indn_p = find(n_new<=1.0e17,1,'first');
        indky_p = find(ky==0,1,'last');
        indn_b = find(n_new<=1.0e17,1,'last');
        indky_b = find(ky<=20,1,'last');
        indn_tm = find(n_new<=1.0e16,1,'last');
        indky_tm = find(ky==0,1,'last');
        indn_bn = find(n_new<=1.0e18,1,'first');
        indky_bn = find(ky>=20,1,'first');
        
        figure(1)
        set(gcf,'Position',[x0 y0 width height],'Color','w')
        ax1 = subplot(2,2,1);
        levelsrkp11 = linspace(0,30,50);
        contourf(log10(n_new),abs(ky),real(kpara11)',levelsrkp11,'Linecolor','none')
        colormap(magma)
        rc11=colorbar;
        rc11.Ruler.TickLabelFormat = '%1.1f';
        rc11.Ruler.Exponent = 1;
        hold on
%         contour(log10(n_new),abs(ky),real(kpara11)',[0.01665,0.01665],'r','Linewidth',2);
%         plot(log10(n_new(indn_np)),ky(indky_np),'sr','MarkerSize',10,'Linewidth',2,...
%             'MarkerFaceColor','r')
%         plot(log10(n_new(indn_p)),ky(indky_p),'.r','MarkerSize',40)
%         plot(log10(n_new(indn_b)),ky(indky_b),'*r','MarkerSize',10,'Linewidth',2)
%         plot(log10(n_new(indn_tm)),ky(indky_tm),'xr','MarkerSize',15,'Linewidth',3)
%         plot(log10(n_new(indn_bn)),ky(indky_bn),'dr','MarkerSize',10,'MarkerFaceColor','r',...
%             'Linewidth',3)
        set(gca,'xtick',[])
%         title('Re[k_{|| 11}] (m^{-1})')
        ylabel('{\itk_{y}} (m^{-1})')
        set(gca,'Fontsize',30,'FontName','CMU Serif')
        set(gca,'Position',[0.15 0.55 0.3 0.35])
        h1 = get(rc11,'Position');
        text(0.02,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
            'color','white')
        hold off

        ax2 = subplot(2,2,2);
        levelsikp11 = linspace(0,2.5,50);
        contourf(log10(n_new),abs(ky),imag(kpara11)',levelsikp11,'Linecolor','none')
        ic11=colorbar();
        colormap(magma)
        ic11.Ruler.TickLabelFormat = '%1.1f';
        hold on
%         contour(log10(n_new),abs(ky),imag(kpara11)',[0.0177,0.0177],'r','Linewidth',2);
        set(gca,'xtick',[])
        set(gca,'ytick',[])
%         title('Im[k_{|| 11}] (m^{-1})')
%         ylabel('k_{y} (m^{-1})')
        set(gca,'Fontsize',30,'FontName','CMU Serif')
        set(gca,'Position',[0.55 0.55 0.3 0.35])
%         set(gca,'colorscale','log')
        h2 = get(ic11,'Position');
        text(0.02,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
            'color','white')
        hold off

        ax3 = subplot(2,2,3);
        levelsrkp21 = linspace(0,0.25,50);
        contourf(log10(n_new),abs(ky),real(kpara21)',levelsrkp21,'Linecolor','none')
%         hold on
%         contourf(log10(n_new),abs(ky),real(kpara21)',[0 0],'Linecolor','k')
        rc21=colorbar();
        rc21.Ruler.Exponent = -1;
%         caxis([0 1])
        rc21.Ruler.TickLabelFormat = '%1.1f';
        colormap(magma)
        hold on
        contour(log10(n_new),abs(ky),real(kpara21)',[0.0027,0.0027],'-.w','Linewidth',3);
%         title('Re[k_{|| 21}] (m^{-1})')
        ylabel('{\itk_{y}} (m^{-1})'),
        xlabel('log_{10}n')
        text(0.02,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
            'color','white')
        set(gca,'Fontsize',30,'FontName','CMU Serif')
%         set(gca,'colorscale','log')
%         h = get(gca,'colorbar');
        set(rc21,'Position',[h1(1) 0.12 h1(3) h1(4)])
        set(gca,'Position',[0.15 0.12 0.3 0.35])

        ax4 = subplot(2,2,4);
        levelsikp21 = linspace(0,35,50);
%         levelsikp21 = logspace(-2,2.1,50);
%         for kk = 1:length(ky)
%             for nn = 1:length(n_new)
%                 if imag(kpara21(nn,kk))<=1.0e-2
%                     kpara21(nn,kk) = real(kpara21(nn,kk)) + 1.0e-2i;
%                 end
%             end
%         end
        contourf(log10(n_new),abs(ky),imag(kpara21)',levelsikp21,'Linecolor','none')
        ic21=colorbar;
        colormap(magma)
        hold on
        ic21.Ruler.TickLabelFormat = '%1.1f';
        ic21.Ruler.Exponent = 1;
        contour(log10(n_new),abs(ky),imag(kpara21)',[0.0013,0.0013],'-.w','Linewidth',3);
%         title('Im[k_{|| 21}] (m^{-1})')
        set(gca,'ytick',[])
%         ylabel('k_{y} (m^{-1})')
        xlabel('log_{10}n')
        set(gca,'Fontsize',30,'FontName','CMU Serif')
        set(gca,'Position',[0.55 0.12 0.3 0.35])
%         set(gca,'colorscale','log')
        set(ic21,'Position',[h2(1) 0.12 h2(3) h2(4)])
        text(0.02,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
            'Color','white')
        hold off
        
%         export_fig('/Volumes/DATA/LAPD/matlab/wave_projection/dispersion_contour_kx10.png',...
%             '-r300')
%         export_fig('/Volumes/DATA/thesis/RFT/figs/dispersion_nonpropcase.png',...
%             '-r300')
        
        figure(3)
        ax1 = subplot(2,2,1);
        levelsrkp11 = logspace(0,2,50);
        contourf(log10(n_new),ky,real(kpara11)','Linecolor','none')
        rc11=colorbar;
        set(gca,'colorscale','log')
%         colormap(flipud(parula))
        set(gca,'xtick',[])
%         title('Re[k_{|| 11}]')
        ylabel('k_{\perp} (m^{-1})')
        set(gca,'Fontsize',20)

        ax2 = subplot(2,2,2);
%         contourf(log10(n_new),ky,real(kpara12)','Linecolor','none')
        contourf(log10(n_new),ky,imag(kpara11)','Linecolor','none')
        rc11 = colorbar;
        set(gca,'colorscale','log')
%         colormap(ax2,flipud(parula))
        set(gca,'xtick',[])
%         title('Re[k_{|| 12}]')
        title('Im[k_{|| 11}]')
        ylabel('k_{\perp} (m^{-1})')
%         xlabel('log_{10}n')
        set(gca,'Fontsize',20)

        ax3 = subplot(2,2,3);
        contourf(log10(n_new),ky,real(kpara21)','Linecolor','none')
        rc11 = colorbar;
%         colormap(flipud(parula))
        set(gca,'colorscale','log')
        title('Re[k_{|| 21}]')
        ylabel('k_{\perp} (m^{-1})')
        xlabel('log_{10}n')
        set(gca,'Fontsize',20)

        ax4 = subplot(2,2,4);
%         contourf(log10(n_new),ky,real(kpara22)','Linecolor','none')
        contourf(log10(n_new),ky,imag(kpara21)','Linecolor','none')
        rc11 = colorbar;
        set(gca,'colorscale','log')
%         colormap(ax4,flipud(parula))
%         title('Re[k_{|| 22}]')
        title('Im[k_{|| 21}]')
        ylabel('k_{\perp} (m^{-1})')
        xlabel('log_{10}n')
        set(gca,'Fontsize',20)
        
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

% figure(9)
% plot(log10(n_new),real(yp11),'.k','Markersize',20)
% 
% hold on
% 
% plot(log10(n_new),imag(yp11),'.r','Markersize',20)
% plot(log10(n_new),real(yp12),'dk','MarkerSize',7)
% plot(log10(n_new),imag(yp12),'dr','MarkerSize',7)
% plot(log10(Nmax)*ones(1,npts),log10(linspace(0.01,100,npts)),'b','Linewidth',3)
% plot(log10(1.0e17)*ones(1,npts),log10(linspace(0.01,100,npts)),'b','Linewidth',3)
% if perp
%     legend('Re[k_{\perp1}]', 'Im[k_{\perp1}]', 'Re[k_{\perp1}]', 'Im[k_{\perp1}]')
%     ylabel('log_{10}|k_{\perp1}|','Fontsize',16)
% elseif para
%     legend('Re[k_{||1}]', 'Im[k_{||1}]', 'Re[k_{||1}]', 'Im[k_{||1}]')
%     ylabel('log_{10}|k_{||1}|','Fontsize',16)
% end
% xlabel('log_{10}|n|','Fontsize',16)
% % vline(log10(Nmin),'--b')
% xlim([log10(min(n_new)),log10(max(n_new))])
% set(gca,'XDir','reverse');
% 
% hold off
% 
% figure(10)
% plot(log10(n_new),real(yp21),'.k','Markersize',20)
% 
% hold on
% 
% plot(log10(n_new),imag(yp21),'.r','Markersize',20)
% plot(log10(n_new),real(yp22),'dk','MarkerSize',7)
% plot(log10(n_new),imag(yp22),'dr','MarkerSize',7)
% if perp
%     legend('Re[k_{\perp2}]', 'Im[k_{\perp2}]', 'Re[k_{\perp2}]', 'Im[k_{\perp2}]')
%     ylabel('log_{10}|k_{\perp2}|','Fontsize',16)
% elseif para
%     legend('Re[k_{||2}]', 'Im[k_{||2}]', 'Re[k_{||2}]', 'Im[k_{||2}]')
%     ylabel('log_{10}|k_{||2}|','Fontsize',16)
% end
% xlabel('log_{10}|n|','Fontsize',16)
% % vline(log10(N0(imme)),'--k') 
% xlim([log10(min(n_new)),log10(max(n_new))])
% set(gca,'XDir','reverse');
% 
% hold off
% 
% 
