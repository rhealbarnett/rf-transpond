
lapd_params;
[om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,npts,1);

plots = 1;

dispersion;

indx1 = find(source(npts/2:npts)==0,1,'first') + npts/2;
indx2 = npts - floor(0.35*npts);
indx = round(linspace(indx1,indx2,indx2-indx1));
NP = length(indx);

kz_spec_density = zeros(npts,npts);

count = 1;

for ii=1
    
    density = n_new(1,ii)*ones(1,npts);
    
    [om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,density,m_s,om,eps0,npts,1);
    [A,rf_e,rf_ex,rf_ey,rf_ez,diss_pow] = wave_sol(zax,ky,kx,k0,...
    om,mu0,cpdt,sig,source,0,1,1);
                                      
    [kz_spec, k_ax, phase] = fft_kz(dx,npts,rf_ex,rf_ey,rf_ez,plots);
%     [kz_spec, k_ax, phase] = fft_kz(dx,NP,rf_ex(indx),rf_ey(indx),rf_ez(indx),plots);
    
    kz_spec_density(count,:) = kz_spec(:,3);
    
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
plot(log10(n_new), real(kp11),'-.r','Linewidth',2)
c = colorbar;
colormap(magma)
caxis([0 1.0e-4])
ylim([0 50])
yticks(linspace(0,50,11))
set(gca,'colorscale','linear')
ylabel('{\it k_z} (m^{-1})')
xlabel('log_{10}({\itn} (m^{-3}))')
ylabel(c,'|FFT[{\it E_z} (Vm^{-1})]|','Fontsize',30)

% export_fig('/Volumes/DATA/thesis/RFT/figs/kvsn_dispersioncont_linearscale.png',...
%     '-r300')

%%

indx = find(n_new<=(10^20),1,'last');

figure(2)
plot(k_ax,kz_spec_density(indx,1:npts/2))