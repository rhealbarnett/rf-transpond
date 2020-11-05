
wave_verification;

[om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,npts,{0,''});
kz_dispersion = dispersion(npts,s_arr,d_arr,p_arr,om,n_refrac,n_new,1,0,ky);

expkz = real(kz_dispersion(1,:));

kz_spec_density = zeros(npts,npts);

count = 1;

for ii=1:npts

    density = n_new(1,ii)*ones(1,npts);

    [om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,density,m_s,om,eps0,npts,{1,damp_len});
    [A,rf_e,rf_ex,rf_ey,rf_ez] = wave_sol(zax,ky,kx,k0,...
    om,mu0,cpdt,source,0,1,1,0);

    [kz_spec, k_ax, phase, dk] = fft_kz(dz,npts,rf_ex,rf_ey,rf_ez,0);

    kz_spec_density(count,:) = kz_spec(:,3);

    ind_kz = find(kz_spec(1:npts/2,2)==max(kz_spec(1:npts/2,2)));
    actkz(1,count) = k_ax(ind_kz);

    count = count + 1;

end

%%

x0 = 0;
y0 = 0;
width = 1000;
height = 500;

indk = find(k_ax<=50);

% for ii=1:npts
%     for jj=1:npts
% 
%         if kz_spec_density(ii,jj)<=1.0e-6
% 
%             kz_spec_density(ii,jj) = 1.0e-6;
% 
%         end
%     end
% end


figure(1)
set(gcf,'Position',[x0 y0 width height],'Color','w')
% levels = logspace(-6,-4,50);
levels = linspace(0,1.e-4,50);
contourf(log10(n_new),k_ax(indk),(kz_spec_density(:,indk))',levels,'Linecolor','none')
hold on
plot(log10(n_new), real(kz_dispersion(1,:)),'-.r','Linewidth',2)
c = colorbar;
colormap(magma)
caxis([0 1.0e-4])
ylim([0 50])
yticks(linspace(0,50,11))
set(gca,'colorscale','linear')
ylabel('{\it k_z} (m^{-1})')
xlabel('log_{10}({\itn} (m^{-3}))')
ylabel(c,'|FFT[{\it E_z} (Vm^{-1})]|','Fontsize',30)

% export_fig('/Volumes/DATA/thesis/figs/gaussian_ne17_epara_fft.png',...
%     '-r300')

%%

indx = find(n_new<=(10^20),1,'last');

figure(2)
plot(k_ax,kz_spec_density(indx,1:npts/2))