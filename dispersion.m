%-----------------------------------------%
% Find kz from quartic, Brambilla sec 18.2%
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

function kz_dispersion = dispersion(npts,s_arr,d_arr,p_arr,om,n_refrac,n_new,para,plots,ky)

    const = constants();
    c0 = const.c0;

    r_arr = s_arr + d_arr;
    l_arr = s_arr - d_arr;
    
    kz_dispersion = zeros(4,npts,length(n_refrac));
    kpara11 = zeros(npts,length(n_refrac));
    kpara12 = zeros(npts,length(n_refrac));
    kpara21 = zeros(npts,length(n_refrac));
    kpara22 = zeros(npts,length(n_refrac));
    
    

    if ~para

        syms kperp
        
        kz = n_refrac*om./c0;

        for ii=1:length(n_refrac)
            a1 = s_arr;
            b1 = r_arr.*l_arr + p_arr.*s_arr - n_refrac(1,ii)^2*(p_arr + s_arr);
            c1 = p_arr.*((n_refrac(1,ii)^2 - r_arr).*(n_refrac(1,ii)^2 - l_arr));

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
%             ns_s = -(n_refrac(1,ii)^2 - s_arr).*(p_arr./s_arr);
%             ns_f = -((n_refrac(1,ii)^2 - r_arr).*(n_refrac(1,ii)^2 - l_arr))./(n_refrac(1,ii)^2 - s_arr);
%             n_s1 = sqrt(ns_s); n_s2 = -sqrt(ns_s); 
%             n_f1 = sqrt(ns_f); n_f2 = -sqrt(ns_f);
%     
%             k_s1 = n_s1*om/c0; k_s2 = n_s2*om/c0; 
%             k_f1 = n_f1*om/c0; k_f2 = n_f2*om/c0;
%             
%             kpara11(:,ii) = k_f1;
%             kpara12(:,ii) = k_f2;
%             kpara21(:,ii) = k_s1;
%             kpara22(:,ii) = k_s2;
    % 
    %         ks1 = sign(k_s1).*log10(abs(k_s1));
    %         ks2 = sign(k_s2).*log10(abs(k_s2));
    %         kf1 = sign(k_f1).*log10(abs(k_f1));
    %         kf2 = sign(k_f2).*log10(abs(k_f2));
        end
        
        kz_dispersion(1,:,:) = kpara11;
        kz_dispersion(2,:,:) = kpara12;
        kz_dispersion(3,:,:) = kpara21;
        kz_dispersion(4,:,:) = kpara22;

        nsize = size(n_new);
        
        k_para = n_refrac*(om/c0);

        if nsize(1)~=1
            n_new = n_new(1,:);
        else
        end

        if length(n_refrac)==1 && plots
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

        elseif length(n_refrac)~=1 && plots
            
            x0 = 0;
            y0 = 0;
            width = 1000;
            height = 800;
            
            
            indn_sp = find(n_new<=1.0e17,1,'last');
            indky_sp = find(kz==3,1,'last');
            indn_spp = find(n_new<=1.0e17,1,'last');
            indky_spp = find(kz<=10,1,'last');
            indn_fp = find(n_new<=1.0e19,1,'last');
            indky_fp = find(kz==3,1,'last');
            indn_fpp = find(n_new<=1.0e19,1,'last');
            indky_fpp = find(kz<=10,1,'last');
            

            figure(1)
            set(gcf,'Position',[x0 y0 width height],'Color','w')
            ax1 = subplot(2,2,1);
            levelsrkp11 = linspace(0,20,512);
            contourf(log10(n_new),abs(k_para),real(kpara11)',levelsrkp11,'Linecolor','none')
            colormap(magma)
%             rc11=colorbar;
%             rc11.Ruler.TickLabelFormat = '%1.1f';
%             rc11.Ruler.Exponent = 1;
%             caxis([0 20])
%             set(rc11,'YTick',[1,1e1])
            hold on
            set(gca,'xtick',[])
    %         title('Re[k_{|| 11}] (m^{-1})')
            ylabel('$k_{z}$ (m$^{-1}$)','interpreter','latex')
            set(gca,'Fontsize',30,'FontName','CMU Serif')
            set(gca,'Position',[0.15 0.55 0.3 0.35])
%             h1 = get(rc11,'Position');
            text(0.02,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','white')
%             set(gca,'colorscale','log')
%             plot(log10(n_new(indn_fp)),kz(indky_fp),'xk','MarkerSize',15,'Linewidth',3)
%             plot(log10(n_new(indn_fpp)),kz(indky_fpp),'xw','MarkerSize',15,'Linewidth',3)
            hold off

            ax2 = subplot(2,2,2);
            levelsikp11 = linspace(0,20,512);
            contourf(log10(n_new),abs(k_para),imag(kpara11)',levelsikp11,'Linecolor','none')
            ic11=colorbar();
            colormap(magma)
%             ic11.Ruler.TickLabelFormat = '%1.1f';
%             ic11.Ruler.Exponent = 1;
            caxis([0 20])
            set(ic11,'YTick',[1,1e1])
            ylabel(ic11,'$k_x$ (m$^{-1}$)','interpreter','latex')
            hold on
    %         contour(log10(n_new),abs(ky),imag(kpara11)',[0.0177,0.0177],'r','Linewidth',2);
            set(gca,'xtick',[])
            set(gca,'ytick',[])
    %         title('Im[k_{|| 11}] (m^{-1})')
    %         ylabel('k_{y} (m^{-1})')
            set(gca,'Fontsize',30,'FontName','CMU Serif')
            set(gca,'Position',[0.55 0.55 0.3 0.35])
            set(gca,'colorscale','log')
            h2 = get(ic11,'Position');
            text(0.02,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','white')
            hold off

            ax3 = subplot(2,2,3);
            levelsrkp21 = linspace(0,1.0e5,2048);
            contourf(log10(n_new),abs(k_para),real(kpara21)',levelsrkp21,'Linecolor','none')
    %         hold on
    %         contourf(log10(n_new),abs(ky),real(kpara21)',[0 0],'Linecolor','k')
            rc21=colorbar();
            caxis([1 10000])
            set(rc21,'YTick',[1e1,1e2,1e3,1e4,1e5])
%             rc21.Ruler.Exponent = 1;
    %         caxis([0 1])
%             rc21.Ruler.TickLabelFormat = '%1.1f';
%             caxis([0 1.0e5])
            colormap(magma)
            hold on
%             contour(log10(n_new),abs(ky),real(kpara21)',[0.0027,0.0027],'-.w','Linewidth',3);
    %         title('Re[k_{|| 21}] (m^{-1})')
            ylabel('$k_{z}$ (m$^{-1}$)','interpreter','latex'),
            xlabel('log$_{10}n$','interpreter','latex')
            text(0.02,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','white')
            set(gca,'Fontsize',30,'FontName','CMU Serif')
            set(gca,'colorscale','log')
    %         h = get(gca,'colorbar');
%             set(rc21,'Position',[h1(1) 0.12 h1(3) h1(4)])
            set(gca,'Position',[0.15 0.12 0.3 0.35])
%             plot(log10(n_new(indn_sp)),kz(indky_sp),'xk','MarkerSize',15,'Linewidth',3)
%             plot(log10(n_new(indn_spp)),kz(indky_spp),'xw','MarkerSize',15,'Linewidth',3)
            hold off

            ax4 = subplot(2,2,4);
            levelsikp21 = linspace(0,1.0e5,2048);
    %         levelsikp21 = logspace(-2,2.1,50);
    %         for kk = 1:length(ky)
    %             for nn = 1:length(n_new)
    %                 if imag(kpara21(nn,kk))<=1.0e-2
    %                     kpara21(nn,kk) = real(kpara21(nn,kk)) + 1.0e-2i;
    %                 end
    %             end
    %         end
            contourf(log10(n_new),abs(k_para),imag(kpara21)',levelsikp21,'Linecolor','none')
            ic21=colorbar;
            colormap(magma)
            caxis([1 10000])
            set(ic21,'YTick',[1e1,1e2,1e3,1e4,1e5])
            ylabel(ic21,'$k_x$ (m$^{-1}$)','interpreter','latex')
            hold on
%             ic21.Ruler.TickLabelFormat = '%1.1f';
%             ic21.Ruler.Exponent = 1;
%             caxis([0 1.0e5])
%             contour(log10(n_new),abs(ky),imag(kpara21)',[0.0013,0.0013],'-.w','Linewidth',3);
    %         title('Im[k_{|| 21}] (m^{-1})')
            set(gca,'ytick',[])
    %         ylabel('k_{y} (m^{-1})')
            xlabel('log$_{10}n$','interpreter','latex')
            set(gca,'Fontsize',30,'FontName','CMU Serif')
            set(gca,'Position',[0.55 0.12 0.3 0.35])
            set(gca,'colorscale','log')
            set(ic21,'Position',[h2(1) 0.12 h2(3) h2(4)])
            text(0.02,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'Color','white')
            hold off

        end

    elseif para

        syms kpara
        
        kpara11 = zeros(npts,length(n_refrac));
        kpara12 = zeros(npts,length(n_refrac));
        kpara21 = zeros(npts,length(n_refrac));
        kpara22 = zeros(npts,length(n_refrac));

        %--
        % initialise kx roots arrays, ensure they are complex
        for ii=1:length(n_refrac)

            a1 = p_arr;
            b1 = (2.0*p_arr.*s_arr - n_refrac(1,ii)^2*(p_arr + s_arr));
            c1 = (n_refrac(1,ii)^2 - p_arr).*(s_arr*n_refrac(1,ii)^2 - r_arr.*l_arr);

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

        end
        
        kz_dispersion(1,:,:) = kpara11;
        kz_dispersion(2,:,:) = kpara12;
        kz_dispersion(3,:,:) = kpara21;
        kz_dispersion(4,:,:) = kpara22;

        if length(n_refrac)==1 && plots

            x0 = 0;
            y0 = 0;
            width = 1000;
            height = 500;

            figure(2)

            set(gcf,'Position',[x0 y0 width height],'Color','w')
            subplot(2,1,1)
            semilogx(n_new, real(kp11),'.k')
            hold on
            semilogx(n_new, real(kp12),'.k')
%             semilogx(n_new, real(kp21),'.k')
%             hold on
%             semilogx(n_new, real(kp22),'.m')
            semilogx(n_new, imag(kp11),'.r')
            semilogx(n_new, imag(kp12),'.r')
            semilogx((1.0e17)*ones(1,npts),(linspace(-30,30,npts)),...
                'b--','Linewidth',1.5)
%             semilogx((1.0e18)*ones(1,npts),(linspace(min(imag(kp22)),max(imag(kp21)),npts)),...
%                 'b--','Linewidth',1.5)
            ylim([-30 30])
%             xlabel('$n$ (m$^{-3}$)','Interpreter','latex')
            xticks([])
            ylabel('$k_{z}$ (m$^{-1}$)','Interpreter','latex')
            set(gca,'Fontsize',25)

            subplot(2,1,2)
            semilogx(n_new, real(kp21),'.k')
            semilogx(n_new, real(kp22),'.k')
%             semilogx(n_new, imag(kp11),'.r')
            hold on
%             semilogx(n_new, imag(kp12),'.b')
            semilogx(n_new, imag(kp21),'.r')
            semilogx(n_new, imag(kp22),'.r')
            semilogx((1.0e17)*ones(1,npts),(linspace(-30,30,npts)),...
                '--b','Linewidth',1.5)
%             semilogx((1.0e18)*ones(1,npts),(linspace(-60,60,npts)),...
%                 '--b','Linewidth',1.5)
            ylim([-30 30])
            xlabel('$n$ (m$^{-3}$)','Interpreter','latex')
            ylabel('$k_{z}$ (m$^{-1}$)','Interpreter','latex')
            
%             yticks([])

            set(gca,'Fontsize',25)

    %         export_fig('/Volumes/DATA/thesis/RFT/figs/kvsn_dispersion_lines.png',...
    %             '-r300')

        elseif length(n_refrac)~=1 && plots

            width = 1000;
            height = 800;
            x0 = 0;
            y0 = 0;

%             indn_np = find(n_new<=1.0e17,1,'first');
%             indky_np = find(ky<=30,1,'last');
            indn_p = find(n_new<=1.0e17,1,'last');
            indky_p = find(ky>=0.0,1,'first');
            indn_pp = find(n_new<=1.0e19,1,'last');
            indky_pp = find(ky<=16.8,1,'last');
            indn_b = find(n_new<=2.0e17,1,'last');
            indky_b = find(ky<=165,1,'last');
            indn_bb = find(n_new<=1.0e17,1,'last');
            indky_bb = find(ky<=580,1,'last');
%             indn_tm = find(n_new<=1.0e16,1,'last');
%             indky_tm = find(ky==0,1,'last');
%             indn_bn = find(n_new<=1.0e18,1,'first');
%             indky_bn = find(ky>=20,1,'first');

            figure(1)
            set(gcf,'Position',[x0 y0 width height],'Color','w')
            ax1 = subplot(2,2,1);
            levelsrkp11 = linspace(0,20,512);
            contourf(log10(n_new),(ky),real(kpara11)',levelsrkp11,'Linecolor','none')
            colormap(magma)
%             rc11=colorbar;
%             rc11.Ruler.TickLabelFormat = '%1.1f';
%             rc11.Ruler.Exponent = 1;
%             caxis([0 25])
            hold on
    %         contour(log10(n_new),abs(ky),real(kpara11)',[0.01665,0.01665],'r','Linewidth',2);
    %         plot(log10(n_new(indn_np)),ky(indky_np),'sr','MarkerSize',10,'Linewidth',2,...
    %             'MarkerFaceColor','r')
%             plot(log10(n_new(indn_p)),ky(indky_p),'xk','MarkerSize',15,'linewidth',8)
%             plot(log10(n_new(indn_pp)),ky(indky_pp),'xw','MarkerSize',15,'linewidth',8)
%             plot(log10(n_new(indn_b)),(ky(indky_b)),'xw','MarkerSize',15,'linewidth',15)
%             plot(log10(n_new(indn_bb)),ky(indky_bb),'xw','MarkerSize',15,'linewidth',8)
            set(gca,'xtick',[])
    %         title('Re[k_{|| 11}] (m^{-1})')
            ylabel('$k_{y}$ (m$^{-1}$)','interpreter','latex')
            set(gca,'Fontsize',30,'FontName','CMU Serif')
            set(gca,'Position',[0.15 0.55 0.3 0.35])
%             h1 = get(rc11,'Position');
            text(0.02,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','white')
            hold off

            ax2 = subplot(2,2,2);
            levelsikp11 = linspace(0,20,512);
            contourf(log10(n_new),(ky),imag(kpara11)',levelsikp11,'Linecolor','none')
            ic11=colorbar();
            ylabel(ic11,'$k_z$ (m$^{-1}$)','interpreter','latex')
            colormap(magma)
%             ic11.Ruler.TickLabelFormat = '%1.1f';
%             ic11.Ruler.Exponent = 1;
            caxis([0 20])
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
            levelsrkp21 = linspace(0,30,512);
            colormap(magma)
%             if sum(real(kpara21(:)))==0
%                 kpara21(1,1) = 1.0 + imag(kpara21(1,1));
%             end
            contourf(log10(n_new),(ky),real(kpara21)',levelsrkp21,'Linecolor','none')
    %         hold on
    %         contourf(log10(n_new),abs(ky),real(kpara21)',[0 0],'Linecolor','k')
%             rc21=colorbar();
%             rc21.Ruler.Exponent = 1;
    %         caxis([0 1])
%             rc21.Ruler.TickLabelFormat = '%1.1f';
%             caxis([0 35])
            
            hold on
%             contour(log10(n_new),abs(ky),real(kpara21)',[0.0027,0.0027],'-.w','Linewidth',3);
    %         title('Re[k_{|| 21}] (m^{-1})')
            ylabel('$k_{y}$ (m$^{-1}$)','interpreter','latex'),
            xlabel('log$_{10}n$','interpreter','latex')
            text(0.02,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','white')
            set(gca,'Fontsize',30,'FontName','CMU Serif')
    %         set(gca,'colorscale','log')
    %         h = get(gca,'colorbar');
%             set(rc21,'Position',[h1(1) 0.12 h1(3) h1(4)])
            set(gca,'Position',[0.15 0.12 0.3 0.35])
            
                

            ax4 = subplot(2,2,4);
            levelsikp21 = linspace(0,30,512);
    %         levelsikp21 = logspace(-2,2.1,50);
    %         for kk = 1:length(ky)
    %             for nn = 1:length(n_new)
    %                 if imag(kpara21(nn,kk))<=1.0e-2
    %                     kpara21(nn,kk) = real(kpara21(nn,kk)) + 1.0e-2i;
    %                 end
    %             end
    %         end
            contourf(log10(n_new),(ky),imag(kpara21)',levelsikp21,'Linecolor','none')
            ic21=colorbar;
            colormap(magma)
            hold on
%             ic21.Ruler.TickLabelFormat = '%1.1f';
%             ic21.Ruler.Exponent = 1;
            caxis([0 30])
            ylabel(ic21,'$k_z$ (m$^{-1}$)','interpreter','latex')
%             contour(log10(n_new),abs(ky),imag(kpara21)',[0.0013,0.0013],'-.w','Linewidth',3);
    %         title('Im[k_{|| 21}] (m^{-1})')
            set(gca,'ytick',[])
    %         ylabel('k_{y} (m^{-1})')
            xlabel('log$_{10}n$','interpreter','latex')
            set(gca,'Fontsize',30,'FontName','CMU Serif')
            set(gca,'Position',[0.55 0.12 0.3 0.35])
    %         set(gca,'colorscale','log')
            set(ic21,'Position',[h2(1) 0.12 h2(3) h2(4)])
            text(0.02,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'Color','black')
            hold off

%             export_fig('/Volumes/DATA/LAPD/matlab/wave_projection/dispersion_contour_kx10.png',...
%                 '-r300')
%             export_fig('/Volumes/DATA/thesis/figs/cmod_dispersion_para_nonprop.png',...
%                 '-r300')

        end

    end
    
end