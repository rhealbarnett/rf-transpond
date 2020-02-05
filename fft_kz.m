
function [kz_spec, k_ax, phase] = fft_kz(dx,npts,rf_ex,rf_ey,rf_ez,plots)

    dk = 1.0/((npts-1)*dx);
    knyq = 1.0/(2.0*dx);
    k_np = (npts/2);
    k_ax = zeros(1,k_np);

    for ii=1:k_np

        k_ax(1,ii) = dk*(ii-1)*2.0*pi;

    end

    fft_Ex = fft(rf_ex)/npts;
    fft_Ey = fft(rf_ey)/npts;
    fft_Ez = fft(rf_ez)/npts;
    
    rf_ex_smoothed = zeros(1,npts);
    rf_ey_smoothed = zeros(1,npts);
    rf_ez_smoothed = zeros(1,npts);

    for ii=2:npts-1
       
        rf_ex_smoothed(1,ii) = 0.25*fft_Ex(1,ii-1) + 0.5*fft_Ex(1,ii) + 0.25*fft_Ex(1,ii+1);
        rf_ey_smoothed(1,ii) = 0.25*fft_Ey(1,ii-1) + 0.5*fft_Ey(1,ii) + 0.25*fft_Ey(1,ii+1);
        rf_ez_smoothed(1,ii) = 0.25*fft_Ez(1,ii-1) + 0.5*fft_Ez(1,ii) + 0.25*fft_Ez(1,ii+1);
        
    end
    
    rf_ex_smoothed(1,1) = fft_Ex(1,1);
    rf_ey_smoothed(1,1) = fft_Ey(1,1);
    rf_ez_smoothed(1,1) = fft_Ez(1,1);
    
    rf_ex_smoothed(1,npts) = fft_Ex(1,npts);
    rf_ey_smoothed(1,npts) = fft_Ey(1,npts);
    rf_ez_smoothed(1,npts) = fft_Ez(1,npts);
    
    kz_spec = [abs(fft_Ex); abs(fft_Ey); abs(fft_Ez)]';
%     kz_spec = [abs(rf_ex_smoothed); abs(rf_ey_smoothed); abs(rf_ez_smoothed)]';

    [mag_ex, idx_ex] = max(abs(fft_Ex(11:npts/2)));
    [mag_ey, idx_ey] = max(abs(fft_Ey(11:npts/2)));
    [mag_ez, idx_ez] = max(abs(fft_Ez(11:npts/2)));

    pex = angle(fft_Ex(idx_ex+10));
    pey = angle(fft_Ey(idx_ey+10));
    pez = angle(fft_Ez(idx_ez+10));

    phase_xy = (pey - pex)*(180/pi);
    phase_xz = (pez - pex)*(180/pi);
    phase_yz = (pey - pez)*(180/pi);

    phase = [phase_xy; phase_xz; phase_yz];

    %%

    if plots
        figure(1)
        subplot(3,1,1)
        plot(k_ax,abs(fft_Ex(1:npts/2)),'k','Linewidth',2)
        xlim([0 100])
        set(gca,'xtick',[])
        ylabel('|FFT[E_x (Vm^{1})]|')

        subplot(3,1,2)
        plot(k_ax,abs(fft_Ey(1:npts/2)),'k','Linewidth',2)
        xlim([0 100])
        set(gca,'xtick',[])
        ylabel('|FFT[E_y (Vm^{1})]|')

        subplot(3,1,3)
        plot(k_ax,abs(fft_Ez(1:npts/2)),'k','Linewidth',2)
        xlim([0 100])
        xlabel('k_z (m^{1})')
        ylabel('|FFT[E_z (Vm^{1})]|')
    end

end