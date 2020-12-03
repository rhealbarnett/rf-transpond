
function [actkz,dk] = kz_spectrum(n_new,q_s,m_s,om,npts,damp_len,dampFac,zax,ky,kx,k0,B0,...
                                source,expkz,plots)
                            
    const = constants();
    
    eps0 = const.eps0;
    mu0 = const.mu0;
        
    dz = zax(2) - zax(1);

    kz_spec_density = zeros(npts,npts);

    count = 1;

    for ii=1:npts

        density = n_new(1,ii)*ones(1,npts);

        [~,~,cpdt,~,~,~,~] = dielec_tens(q_s,B0,density,m_s,om,...
            eps0,npts,{1,damp_len,dampFac});
        [~,~,rf_ex,rf_ey,rf_ez] = rf_wave_sol(zax,ky,kx,k0,...
        om,mu0,cpdt,source,0,1,1);

        [kz_spec, k_ax, ~, dk] = fft_kz(dz,npts,rf_ex,rf_ey,rf_ez,0);

        kz_spec_density(count,:) = kz_spec(:,3);

        ind_kz = find(kz_spec(1:npts/2,2)==max(kz_spec(1:npts/2,2)));
        actkz(1,count) = k_ax(ind_kz);

        count = count + 1;

    end

    if plots
        x0 = 0;
        y0 = 0;
        width = 1000;
        height = 500;

        indk = find(k_ax<=50);

        figure(5)
        set(gcf,'Position',[x0 y0 width height],'Color','w')
        levels = linspace(0,5.e-5,50);
        contourf(log10(n_new),k_ax(indk),(kz_spec_density(:,indk))',levels,'Linecolor','none')
        hold on
        plot(log10(n_new), expkz,'-.r','Linewidth',2)
        c = colorbar;
        colormap(magma)
        caxis([0 5.0e-5])
        ylim([0 50])
        yticks(linspace(0,50,11))
        set(gca,'colorscale','linear')
        ylabel('{\it k_z} (m^{-1})')
        xlabel('log_{10}({\itn} (m^{-3}))')
        ylabel(c,'|FFT[{\it E_z} (Vm^{-1})]|','Fontsize',30)
        
        saveas(gcf,'outputs/wave_verification.png');
        close 5
    end

end
