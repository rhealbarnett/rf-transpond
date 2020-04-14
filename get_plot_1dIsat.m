%% Call in xz plane Isat data file

filename = '/Volumes/DATA/LAPD/RF_April22/07_Mach_p30y_12kV.hdf5';


%% get Mach probe data: check the LAPD_201904_documentation ---------------- %%


channels = ["Mach1", "Mach2", "Mach3", "Mach4","Mach5","Mach6"];


% Mach 1 & 2: x, orientation: if 1>2, flow is -x (toward antenna)
% Mach 3 & 4: y, orientation: if 3>4, flow is -y (vertically down)
% Mach 5 & 6: z, orientation: if 5>6, flow is -z (toward LaB6 source)
% To see an approximation for the density, average the paired probes.
% To see an approximation for the flow, take the difference of the paired
% probes. 

z_line = 1;

% attenuation (V) and resistance (omhs), needed to calibrate the signal.
atten = 3.16;
res = 25;

% time values taken 'safely' within the rf application.
rfmin = 0.7;
rfmax = 1.2;

% time values before the rf appliation.
tmin = 0.05;
tmax = 0.5;

plots = 1;

%%

if z_line
    % xz plane near antenna
    % (-11 <= x <= -6) cm, (-12 <= z <= 12) cm
    
    filename = '/Volumes/DATA/LAPD/RF_April22/11_Mach_p30xyz_12kV.hdf5';
    
    probename = 'Mach';
    
    gcx = -8.5;
    gcz = 0.0;

    dx = 0.5;
    dz = 0.25;

    Nx = 11;
    Nz = 97;

    minx = gcx - ((Nx-1)/2)*dx;
    maxx = gcx + ((Nx-1)/2)*dx;

    minz = gcz - ((Nz-1)/2)*dz;
    maxz = gcz + ((Nz-1)/2)*dz;

    % x = linspace(minx,maxx,Nx);
    x = -11;
    z = linspace(minz,maxz,Nz);


     % get data ---------------------------------------------------------------- %

    [actx_LP, acty_LP, time_LP, data_LP] = get_probe(filename, probename, channels, x, z, dz, 1);

    callib = atten/res;

    for ii= 1:numel(x)
        for jj= 1:numel(z)
            data_LP{ii,jj,1} = atten*data_LP{ii,jj,1} / res; % multiply by attenuated value, divide by resistance
            data_LP{ii,jj,2} = atten*data_LP{ii,jj,2} / res; 
            data_LP{ii,jj,3} = atten*data_LP{ii,jj,3} / res; 
            data_LP{ii,jj,4} = atten*data_LP{ii,jj,4} / res;
            data_LP{ii,jj,5} = atten*data_LP{ii,jj,5} / res;
            data_LP{ii,jj,6} = atten*data_LP{ii,jj,6} / res;
        end
    end

elseif ~z_line
    
    probename = 'XY[3]: 3-axisMach_Suying#1';
    
    gcx = -10.0;
    gcy = 0.0;

    dx = 0.0;
    dy = 1.0;

    Nx = 1;
    Ny = 21;

    minx = gcx - ((Nx-1)/2)*dx;
    maxx = gcx + ((Nx-1)/2)*dx;

    miny = gcy - ((Ny-1)/2)*dy;
    maxy = gcy + ((Ny-1)/2)*dy;

    % x = linspace(minx,maxx,Nx);
    x = -10;
    y = linspace(miny,maxy,Ny);
    
    temp_x = zeros(120832,Ny,12);
    temp_y = zeros(120832,Ny,12);
    temp_z = zeros(120832,Ny,12);
    
    fft_sigx = zeros(120832,Ny,12);
    fft_sigy = zeros(120832,Ny,12);
    fft_sigz = zeros(120832,Ny,12);

    for kk=1:12

        filename = ['/Volumes/DATA/LAPD/RF_April22/07_Mach_p30y_', num2str(kk), 'kV.hdf5'];

        [actx_LP, acty_LP, time_LP, data_LP] = get_probe(filename, probename, channels, x, y, dy/2, 1);

%     [actx_LP, acty_LP, time_LP, data_LP] = get_probe(filename, probename, channels, x, y, dy, 1);
    
        for ii= 1:numel(x)
            for jj= 1:numel(y)
                data_LP{ii,jj,1} = atten*data_LP{ii,jj,1} / res; % multiply by attenuated value, divide by resistance
                data_LP{ii,jj,2} = atten*data_LP{ii,jj,2} / res; 
                data_LP{ii,jj,3} = atten*data_LP{ii,jj,3} / res; 
                data_LP{ii,jj,4} = atten*data_LP{ii,jj,4} / res;
                data_LP{ii,jj,5} = atten*data_LP{ii,jj,5} / res;
                data_LP{ii,jj,6} = atten*data_LP{ii,jj,6} / res;
            end
        end

        for ii= 1:numel(x)
            for jj= 1:numel(y)
                data_x{ii,jj,1} = (data_LP{ii,jj,1} + data_LP{ii,jj,2})./2.0;
                data_y{ii,jj,1} = (data_LP{ii,jj,3} + data_LP{ii,jj,4})./2.0;
                data_z{ii,jj,1} = (data_LP{ii,jj,5} + data_LP{ii,jj,6})./2.0;

                temp_x(:,jj,kk) = data_x{ii,jj,1};
                temp_y(:,jj,kk) = data_y{ii,jj,1};
                temp_z(:,jj,kk) = data_z{ii,jj,1};

            end
        end
        
    end
    
    npts = length(time_LP);
    dt = 4.0000e-08;
    df = 1.0/((npts)*dt);
    Fnyq = 1.0/(2.0*dt);
    freq_np = (npts/2);
    freq_ax = zeros(1,freq_np);

    for ii=1:freq_np
        freq_ax(1,ii) = df*(ii-1);
    end
    
    bins = find(freq_ax > 1e5);
    
    for ii = 1:numel(y)
        for jj = 1:12
    
            fft_sigx(:,ii,jj) = fft(temp_x(:,ii,jj));
            fft_sigy(:,ii,jj) = fft(temp_y(:,ii,jj));
            fft_sigz(:,ii,jj) = fft(temp_z(:,ii,jj));
            
            fft_sigx(int32(bins),ii,jj) = 0.0;
            fft_sigx(npts - int32((bins)),ii,jj) = 0.0;
            inv_sigx(:,ii,jj) = ifft(fft_sigx(:,ii,jj));
            filt_x(:,ii,jj) = real(inv_sigx(:,ii,jj));
            fft_sigy(int32(bins),ii,jj) = 0.0;
            fft_sigy(npts - int32((bins)),ii,jj) = 0.0;
            inv_sigy(:,ii,jj) = ifft(fft_sigy(:,ii,jj));
            filt_y(:,ii,jj) = real(inv_sigy(:,ii,jj));
            fft_sigz(int32(bins),ii,jj) = 0.0;
            fft_sigz(npts - int32((bins)),ii,jj) = 0.0;
            inv_sigz(:,ii,jj) = ifft(fft_sigz(:,ii,jj));
            filt_z(:,ii,jj) = real(inv_sigz(:,ii,jj));
            
        end
    end
    
end


%% Calculate average of each probe for density proxy from Isat


for ii= 1:numel(x)
    for jj= 1:numel(z)
        data_x{ii,jj,1} = (data_LP{ii,jj,1} + data_LP{ii,jj,2})./2.0;
        data_y{ii,jj,1} = (data_LP{ii,jj,3} + data_LP{ii,jj,4})./2.0;
        data_z{ii,jj,1} = (data_LP{ii,jj,5} + data_LP{ii,jj,6})./2.0;
    end
end

temp = data_x{1,floor(Nz/2),1};

x0 = 0;
y0 = 0;
width = 1000;
height = 500;

figure(1)
set(gcf,'Position',[x0 y0 width height],'Color','w')
plot(time_LP(1:10:end),temp(1:10:end),'--k','Linewidth',0.5)
ylabel('I_{\it sat} (A)')
xlabel('Time (ms)')
xlim([min(time_LP) max(time_LP)])
hold on

%% Find Isat > 4.0ms, average, and shift data by this value. 

shift_time = find(time_LP>=(max(time_LP)-1.0));

for ii= 1:numel(x)
    for jj= 1:numel(z)
        tempx = data_x{ii,jj,1};
        tempy = data_y{ii,jj,1};
        tempz = data_z{ii,jj,1};
        data_x{ii,jj,1} = data_x{ii,jj,1} + abs(mean(tempx(shift_time)));
        data_y{ii,jj,1} = data_y{ii,jj,1} + abs(mean(tempy(shift_time)));
        data_z{ii,jj,1} = data_z{ii,jj,1} + abs(mean(tempz(shift_time)));
    end
end

clear tempx tempy tempz

temp = data_x{1,floor(Nz/2),1};

figure(1)
plot(time_LP,temp,'k','Linewidth',1.0)

% export_fig('/Volumes/DATA/thesis/RFT/figs/isat_timeseries.png',...
%     '-r300')

%% Calculate moving mean for the raw Isat data

for ii= 1:numel(x)
    for jj= 1:numel(z)
        mm_x{ii,jj,1} = movmean(data_x{ii,jj,1},121);
        mm_y{ii,jj,1} = movmean(data_y{ii,jj,1},121);
        mm_z{ii,jj,1} = movmean(data_z{ii,jj,1},121);
    end
end

%% Plot raw signal and FFT of raw signal

t = find(time_LP > 0.49 & time_LP < 0.51);
npts = length(time_LP);
tax = time_LP*1.0e-3;
% temp = data_x{1,(Nz+1)/2,1};
plot_temp = mm_x{1,floor(Nz/2),1};

%% Plot moving mean results.

figure(1)
subplot(2,2,1)
plot(tax,temp)
xlabel('Time (s)')
ylabel('Raw Isat (A)')
title('Raw Isat')

subplot(2,2,3)
plot(tax(t),temp(t))
hold on
plot(0.501e-3*ones(1,length(tax)),linspace(-0.02,0.05,length(tax)),'--k')
xlabel('Time (s)')
ylabel('Raw Isat (A)')
title('Raw Isat Zoomed')

subplot(2,2,2)
plot(tax,plot_temp)
xlabel('Time (s)')
ylabel('Filtered Isat (A)')
title('Filtered Isat (moving mean, wind length 121)')

subplot(2,2,4)
plot(tax(t),plot_temp(t))
hold on
plot(0.501e-3*ones(1,length(tax(t))),linspace(0.012,0.018,length(tax(t))),'--k')
xlabel('Time (s)')
ylabel('Filtered Isat (A)')
title('Filtered Isat Zoomed')


%% FFT of raw Isat data

npts = length(time_LP);
dt = 4.0000e-08;
df = 1.0/((npts)*dt);
Fnyq = 1.0/(2.0*dt);
freq_np = (npts/2);
freq_ax = zeros(1,freq_np);

for ii=1:freq_np
    freq_ax(1,ii) = df*(ii-1);
end

data_mean = mean(temp);
fft_sig = (fft(temp - data_mean));

high_freq = find(freq_ax > 1.0e4 & freq_ax < Fnyq);
pow = 2.0*(fft_sig(1:npts/2).*conj(fft_sig(1:npts/2)))*(1.0/npts);
freq_data_ind = find(pow(high_freq)==max(pow(high_freq)));
freq_data = freq_ax(freq_data_ind);
T_data = 1.0/freq_data;

%%

%% Plot first portion of the FFT data

figure(2)
set(gcf,'Position',[x0 y0 width height],'color','w')
% subplot(2,2,1)
% plot(tax,temp)
% xlabel('Time (s)')
% ylabel('Raw Isat (A)')
% title('Raw Isat')

subplot(1,2,1)
plot(freq_ax,pow,'k','Linewidth',2.0)
ax1 = get(gca);
xlabel('Frequency (Hz)')
ylabel('|FFT[I_{\it sat} (A)]|')
xlim([min(freq_ax)-1.0e3*df max(freq_ax)+1.0e2*df])
xticks([0.0, 2.5, 5.0, 7.5, 10.0, 12.5]*1.0e6)
% title('FFT of raw Isat')


hold on

%% Find and zero FFT values in frequency bins above some value. 

bins = find(freq_ax > 1e5);
% fft_sig(int32(bins)) = 0.0;
% fft_sig(npts - int32((bins))) = 0.0;
% inv_sig = ifft(fft_sig);

% filt_x = real(inv_sig);

for ii= 1:numel(x)
    for jj= 1:numel(z)
        tempx = data_x{ii,jj,1};
        tempy = data_y{ii,jj,1};
        tempz = data_z{ii,jj,1};
        
        datax_mean = mean(tempx);
        fft_sigx = (fft(tempx));
        datay_mean = mean(tempy);
        fft_sigy = (fft(tempy));
        dataz_mean = mean(tempz);
        fft_sigz = (fft(tempz));
        
        fft_sigx(int32(bins)) = 0.0;
        fft_sigx(npts - int32((bins))) = 0.0;
        inv_sigx = ifft(fft_sigx);
        filt_x = real(inv_sigx);
        fft_sigy(int32(bins)) = 0.0;
        fft_sigy(npts - int32((bins))) = 0.0;
        inv_sigy = ifft(fft_sigy);
        filt_y = real(inv_sigy);
        fft_sigz(int32(bins)) = 0.0;
        fft_sigz(npts - int32((bins))) = 0.0;
        inv_sigz = ifft(fft_sigz);
        filt_z = real(inv_sigz);
        
        fft_datax{ii,jj,1} = filt_x;
        fft_datay{ii,jj,1} = filt_y;
        fft_dataz{ii,jj,1} = filt_z;
    end
end

clear tempx tempy tempz

%% Plots of raw Isat data and FFTs

zprime = z + 1;
shift_time = find(time_LP<=(max(time_LP)-1.0));
indx = find(freq_ax <= 3.0e6);
indz_ant = find(zprime==0);
indz_edge = find(zprime==13);
temp_ant = data_x{1,indz_ant,1};
temp_edge = data_x{1,indz_edge,1};
% indx = bins(1);
temp = fft_datax{1,floor(Nz/2),1};
rf_zoom = find(time_LP > 0.4 & time_LP < 0.65);
freq = 2.52e6;
TRF = 1.0/freq;
TRF_ax = find(time_LP <= 0.5008+(100*TRF*1.0e3),1,'last');

figure(2)

subplot(1,2,1)
ax1 = get(gca);
hold on
ax2 = axes('Position',[.25 .7 .2 .2]);
box on
plot(freq_ax(indx),pow(indx),'k','Linewidth',1.5);
hold on
ylim([0 18])
plot(bins(end)*ones(1,length(freq_ax)),linspace(0,18,length(freq_ax)),'--r',...
    'Linewidth',1.5)
yticks([])

subplot(1,2,2)
ax1 = plot(time_LP(shift_time),temp(shift_time),'k','Linewidth',1.5);
hold on
plot(0.5008*ones(1,length(temp(shift_time))),...
    linspace(0.0,0.025,length(temp(shift_time))),'--r','Linewidth',1.5,...
    'HandleVisibility','off')
xlabel('Time (ms)')
ylabel('I_{\it sat} (A)')
xlim([min(time_LP) max(time_LP(shift_time))])
ax2 = axes('Position',[.8 .7 .08 .2]);
box on
plot(time_LP(rf_zoom),temp(rf_zoom),'-k','Linewidth',1.5)
hold on
plot(0.5008*ones(1,length(temp(shift_time))),...
    linspace(0.013,0.02,length(temp(shift_time))),'--r','Linewidth',1.5,...
    'HandleVisibility','off')
% plot(time_LP(TRF_ax)*ones(1,length(temp(shift_time))),...
%     linspace(0.013,0.02,length(temp(shift_time))),'--b','Linewidth',1.5,...
%     'HandleVisibility','off')
yticks([])
% title('Filtered Isat')

% subplot(2,2,4)
% plot(tax(t),filt_x(t))
% xlabel('Time (s)')
% ylabel('Filtered Isat (A)')
% hold on
% plot(0.501e-3*ones(1,length(tax(t))),linspace(0.012,0.018,length(tax(t))),'--k')
% title('Filtered Isat Zoomed')
hold off

% export_fig('/Volumes/DATA/thesis/RFT/figs/isat_fft_filtered.png',...
%     '-r300')

x0 = 0;
y0 = 0;
width = 500;
height = 600;
figure(3)
set(gcf,'Position',[x0 y0 width height],'Color','w')
plot(time_LP(shift_time),temp_ant(shift_time),'k','Linewidth',1.5)
hold on
plot(time_LP(shift_time),temp_edge(shift_time),'b','Linewidth',1)
legend('0 cm (antenna)','13 cm')
plot(0.5008*ones(1,length(temp(shift_time))),...
    linspace(0.0,0.025,length(temp(shift_time))),'--r','Linewidth',1.5,...
    'HandleVisibility','off')
xlabel('Time (ms)')
ylabel('I_{\it sat} (A)')
xlim([min(time_LP) max(time_LP(shift_time))])


% export_fig('/Volumes/DATA/LAPD/isat_ant_edge.png',...
%     '-r300')

%% Find indices for each number of RF periods

T_ms = T_data*1.0e3;

tRF = 0.501e-3;
t4 = 100*T_data;
tend = 1.45e-3;

tRF_ax = find(time_LP<=(tRF),1,'last');
t4_ax = find(time_LP<=(tRF+t4),1,'last');

%% Plot 

figure(3)
plot(tax,filt_x)
hold on
plot(tRF*ones(1,length(time_LP)),linspace(-0.01,0.025,length(time_LP)),'--k')
plot((tRF+t4)*ones(1,length(time_LP)),linspace(-0.01,0.025,length(time_LP)),'--r')
plot(tend*ones(1,length(time_LP)),linspace(-0.01,0.025,length(time_LP)),'--b')
legend('I_{sat}','t_{RF}', 't_{RF} +100T', 't_{RF,end}')


%% Linear fit data before RF (approx 0.05ms - 0.4ms)

t_before = find(tax > tmin*1.0e-3 & tax < tmax*1.0e-3);
tind4 = find(tax > tRF+t4 & tax < tend);

coeffs_before = polyfit(tax(t_before),filt_x(t_before),1);
fit_before = polyval(coeffs_before,tax(t_before));

coeffs_t4 = polyfit(tax(tind4),filt_x(tind4),1);
fit_t4 = polyval(coeffs_t4,tax(tind4));

figure(3)
hold on
plot(tax(t_before),fit_before,'Linewidth',3)
plot(tax(tind4),fit_t4,'Linewidth',3)
set(0,'DefaultLegendAutoUpdate','off')
text(tax(t_before(floor(end/2))),fit_before(floor(end/2))+0.1*fit_before(floor(end/2)),...
    ['y = ' num2str(coeffs_before(1)) '*x + ' num2str(coeffs_before(2))])
text(tax(tind4(floor(end/2))),fit_t4(floor(end/2))+0.1*fit_t4(floor(end/2)),...
    ['y = ' num2str(coeffs_t4(1)) '*x + ' num2str(coeffs_t4(2))])
hold off



%% 

filt_x = zeros(Nz,length(tax));

for ii= 1:numel(x)
    for jj=1:numel(z)
        
        temp = data_x{ii,jj,1};
        fft_sig = (fft(temp));
        
        bins = find(freq_ax > 1e5);
        fft_sig(int32(bins)) = 0.0;
        fft_sig(npts - int32((bins))) = 0.0;
        inv_sig = ifft(fft_sig);

        filt_x(jj,:) = real(inv_sig);
            
    end
end

%%

ind = find(z <= 4.0 & z>=-4.0);

for ii=ind(1):4:ind(end)
    
    figure(4)
    plot(tax,filt_x(ii,:),'DisplayName', ['z = ' num2str(z(ii)) ' cm'])
    hold on
    xlabel('Time (s)')
    ylabel('Filtered Isat (A)')
    legend('show')
    
end

figure(4)
hold on
set(gca,'Fontsize',20)
hold off

%%

T_ms = T_data*1.0e3;

tRF = 0.5008;
t4 = 100*T_ms;

tRF_ax = find(time_LP<=(tRF),1,'last');
t4_ax = find(time_LP<=(tRF+t4),1,'last');

t_before = find(time_LP > tmin & time_LP < tmax);
% tind4 = find(time_LP > tRF+t4 & time_LP < tend);
tafter_npts = length(t_before);


for ii=1:Nz
    
    tempz = fft_dataz{1,ii,1};
    tempy = fft_datay{1,ii,1};
    tempx = fft_datax{1,ii,1};
    
    zavg_before = mean(tempz(t_before));
    zavg_after = mean(tempz(t4_ax:(t4_ax+tafter_npts)));
    yavg_before = mean(tempy(t_before));
    yavg_after = mean(tempy(t4_ax:(t4_ax+tafter_npts)));
    xavg_before = mean(tempx(t_before));
    xavg_after = mean(tempx(t4_ax:(t4_ax+tafter_npts)));
    
%     ratio_z(1,ii) = (temp(tRF_ax) - temp(t4_ax))/temp(tRF_ax);
    ratio_z(1,ii) = (zavg_before - zavg_after) / zavg_before;
    ratio_y(1,ii) = (yavg_before - yavg_after) / yavg_before;
    ratio_x(1,ii) = (xavg_before - xavg_after) / xavg_before;
    
end

ratio = (ratio_x + ratio_y + ratio_z)/3.;

figure(5)
set(gcf,'Position',[x0 y0 width height],'color','w')
set(gca,'Fontsize',30)
plot(z,ratio_x,'-.k','Linewidth',1.5)
hold on
plot(z,ratio_y,'-ok','Linewidth',1.5)
plot(z,ratio_z,'-*k','Linewidth',1.5)
legend('R_{I{\itsat,x}}','R_{I{\itsat,y}}','R_{I{\itsat,z}}','location','north west')
xlim([min(z) max(z)])
xlabel('Axial Position (cm)')
ylabel('(I_{\itsat} - I_{{\itsat},RF})/I_{\itsat}')

% export_fig('/Volumes/DATA/thesis/RFT/figs/isat_ratio_line.png',...
%     '-r300')

figure(6)
set(gcf,'Position',[x0 y0 width height],'color','w')
set(gca,'Fontsize',30)
plot(z,ratio,'-*k','Linewidth',1.5)
hold on
% legend('R_{I{\itsat,x}}','R_{I{\itsat,y}}','R_{I{\itsat,z}}','location','north west')
xlim([min(z) max(z)])
xlabel('$z$ (cm)','Interpreter','latex')
ylabel('{\itR}_{I_{\itsat}}')
plot(-4*ones(1,Nz),linspace(0.2,0.4,Nz),...
            '--b','Linewidth',1.5,'HandleVisibility','off')
plot(2*ones(1,Nz),linspace(0.2,0.4,Nz),...
            '--b','Linewidth',1.5,'HandleVisibility','off')
ylim([0.2,0.4])

% export_fig('/Users/rhealbarnett/Documents/Documents/presentations/2020-rfscidac/isat_ratio_avg_line.png',...
%     '-r300')

%% Isat image plots
    
figure(6)
subplot(1,3,1)
imagesc('XData', actx_LP(:,1), 'YData', acty_LP(1,:), 'CData', data_2x');
ylabel('Axial direction (cm)'); xlabel('Radial direction (cm)');
% c=colorbar; colormap(flipud(hot)); 
lim = caxis; caxis([0 0.04])
% xlim = ([minx-0.5*dx maxx+0.5*dx]);
set(gca, 'XDir','reverse')
title('Isat x (A?)')
% set(gca, 'XTickLabel', [])

subplot(1,3,2)
imagesc('XData', actx_LP(:,1), 'YData', acty_LP(1,:), 'CData', data_2y');
xlabel('Radial direction (cm)');
% c=colorbar; colormap(flipud(hot)); 
lim = caxis; caxis([0,0.04])
% xlim = ([minx-0.5*dx maxx+0.5*dx]);
set(gca, 'XDir','reverse')
% set(gca, 'XTickLabel', [])
title('Isat y (A?)')
set(gca, 'YTickLabel', [])

subplot(1,3,3)
imagesc('XData', actx_LP(:,1), 'YData', acty_LP(1,:), 'CData', data_2z');
xlabel('Radial direction (cm)'); 
c=colorbar; colormap(flipud(hot)); 
title('Isat z (A?)')
lim = caxis; caxis([0,0.04])
% xlim = ([minx-0.5*dx maxx+0.5*dx]);
set(gca, 'XDir','reverse')
set(gca, 'YTickLabel', [])


%% Isat line plots over z

figure(7)
plot(z,data_2x(1,:),'.-')

hold on

plot(z,data_T1dx,'.-')
plot(z,data_T2dx,'.-')
plot(z,data_T3dx,'.-')
plot(z,data_T4dx,'.-')

xlim([minz maxz])
xlabel('z Position (cm)')
ylabel('Saturation current (A?)')
title('Line taken at x = -11cm (approx 1cm in front of antenna, located at z=(0\pm3)cm)')
legend('Before RF','25T','50T','75T','100T')

hold off














