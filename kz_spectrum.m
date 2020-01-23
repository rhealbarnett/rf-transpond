
lapd_params;
[om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,npts,0);

plots = 0;

dispersion;

kz_spec_density = zeros(floor(npts/10),npts);

count = 1;

for ii=1:npts
    
    density = n_new(1,ii)*ones(1,npts);
    
    [om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,density,m_s,om,eps0,npts,0);
    [A,rf_e,rf_ex,rf_ey,rf_ez,diss_pow] = wave_sol(zax,ky,kx,k0,...
    om,mu0,cpdt,source,0,1,1);
    
    [kz_spec, k_ax, phase] = fft_kz(dx,npts,rf_ex,rf_ey,rf_ez,plots);
    
    kz_spec_density(count,:) = kz_spec(:,3);
    
    count = count + 1;
    
end

%%

x0 = 0;
y0 = 0;
width = 600;
height = 600;

indk = find(k_ax<=50);

figure(1)
set(gcf,'Position',[x0 y0 width height],'Color','w')
levels = logspace(-4,4,50);
contourf(log10(n_new),k_ax(indk),(kz_spec_density(:,indk))',levels,'Linecolor','none')
hold on
plot(log10(n_new), real(kp11),'k')
c = colorbar;
% colormap(flipud(gray))
caxis([1.0e-4 1.0e4])
set(gca,'colorscale','log')
ylabel('k_z (m^{-1})')
xlabel('log_{10}(n (m^{-3}))')
ylabel(c,'|FFT[E_z (Vm^{-1})]|','Fontsize',20)

% export_fig('/Volumes/DATA/matlab/wave_verification/kvsn_dispersioncont.png',...
%     '-r300')