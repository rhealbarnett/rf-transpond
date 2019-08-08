%% Call in xz plane Isat data file

filename = '/Users/rhealbarnett/Documents/LAPD/RF_April22/11_Mach_p30xyz_12kV.hdf5';


%% get Mach probe data: check the LAPD_201904_documentation ---------------- %%

probename = 'Mach';
channels = ["Mach1", "Mach2", "Mach3", "Mach4","Mach5","Mach6"];


% Mach 1 & 2: x, orientation: if 1>2, flow is -x (toward antenna)
% Mach 3 & 4: y, orientation: if 3>4, flow is -y (vertically down)
% Mach 5 & 6: z, orientation: if 5>6, flow is -z (toward LaB6 source)
% To see an approximation for the density, average the paired probes.
% To see an approximation for the flow, take the difference of the paired
% probes. 

% xz plane near antenna
% (-11 <= x <= -6) cm, (-12 <= z <= 12) cm
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

% attenuation (V) and resistance (omhs), needed to calibrate the signal.
atten = 3.16;
res = 25;

% time values taken 'safely' within the rf application.
rfmin = 0.7;
rfmax = 1.2;

% time values before the rf appliation.
tmin = 0.05;
tmax = 0.4;

plots = 1;


%% get data ---------------------------------------------------------------- %

[actx_LP, acty_LP, time_LP, data_LP] = get_probe(filename, probename, channels, x, z, 1);

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


%% Calculate average of each probe for density proxy from Isat


for ii= 1:numel(x)
    for jj= 1:numel(z)
        data_x{ii,jj,1} = (data_LP{ii,jj,1} + data_LP{ii,jj,2})./2.0;
        data_y{ii,jj,1} = (data_LP{ii,jj,3} + data_LP{ii,jj,4})./2.0;
        data_z{ii,jj,1} = (data_LP{ii,jj,5} + data_LP{ii,jj,6})./2.0;
    end
end

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
temp = data_x{1,(Nz+1)/2,1};
plot_temp = mm_x{1,(Nz+1)/2,1};

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

dt = 4.0000e-08;
df = 1.0/((npts)*dt);
Fnyq = 1.0/(2.0*dt);
freq_np = (npts/2);
freq_ax = zeros(1,freq_np);

for ii=1:freq_np
    freq_ax(1,ii) = df*(ii-1);
end

data_mean = mean(temp);
fft_sig = (fft(temp));

%% Plot first portion of the FFT data

figure(2)
subplot(2,2,1)
plot(tax,temp)
xlabel('Time (s)')
ylabel('Raw Isat (A)')
title('Raw Isat')

subplot(2,2,3)
plot(freq_ax,2.0*(fft_sig(1:npts/2).*conj(fft_sig(1:npts/2)))*(1.0/npts))
xlabel('Frequency (Hz)')
ylabel('2*abs[FFT(Raw Isat)]')
title('FFT of raw Isat')

hold on

%% Find and zero FFT values in frequency bins above some value. 

bins = find(freq_ax > 1e5);
fft_sig(int32(bins)) = 0.0;
fft_sig(npts - int32((bins))) = 0.0;
inv_sig = ifft(fft_sig);

filt_x = real(inv_sig);

%% Plots of raw Isat data and FFTs

figure(2)

subplot(2,2,2)
plot(tax,filt_x)
xlabel('Time (s)')
ylabel('Filtered Isat (A)')
title('Filtered Isat')

subplot(2,2,4)
plot(tax(t),filt_x(t))
xlabel('Time (s)')
ylabel('Filtered Isat (A)')
hold on
plot(0.501e-3*ones(1,length(tax(t))),linspace(0.012,0.018,length(tax(t))),'--k')
title('Filtered Isat Zoomed')
hold off


%% Find indices for each number of RF periods

f = 2.38e6;
T = 1.0/f;
T_ms = T*1.0e3;

TRF = 0.501;
T1 = 25*T_ms;
T2 = 50*T_ms;
T3 = 75*T_ms;
T4 = 100*T_ms;

TRF_ax = find(time_LP<=(TRF),1,'last');
T1_ax = find(time_LP<=(TRF+T1),1,'last');
T2_ax = find(time_LP<=(TRF+T2),1,'last');
T3_ax = find(time_LP<=(TRF+T3),1,'last');
T4_ax = find(time_LP<=(TRF+T4),1,'last');

%% Plot 

figure(3)
plot(time_LP,filt_x)
hold on
plot(TRF*ones(1,length(time_LP)),linspace(-0.01,0.025,length(time_LP)),'--k')
plot((TRF+T1)*ones(1,length(time_LP)),linspace(-0.01,0.025,length(time_LP)),'--r')
plot((TRF+T2)*ones(1,length(time_LP)),linspace(-0.01,0.025,length(time_LP)),'--b')
plot((TRF+T3)*ones(1,length(time_LP)),linspace(-0.01,0.025,length(time_LP)),'--c')
plot((TRF+T4)*ones(1,length(time_LP)),linspace(-0.01,0.025,length(time_LP)),'--m')
legend('I_{sat}','T_{RF}','T_{RF} +25T', 'T_{RF} +50T', 'T_{RF} +75T', 'T_{RF} +100T')
% plot((0.501-T1)*ones(1,length(time_LP)),linspace(-0.005,0.02,length(time_LP)),'--r')
% plot((0.501-T2)*ones(1,length(time_LP)),linspace(-0.005,0.02,length(time_LP)),'--b')
% plot((0.501-T3)*ones(1,length(time_LP)),linspace(-0.005,0.02,length(time_LP)),'--c')
% plot((0.501-T4)*ones(1,length(time_LP)),linspace(-0.005,0.02,length(time_LP)),'--m')
hold off


%% 

for ii= 1:numel(x)
    for jj=1:numel(z)
        
        t_Langmuir = find(time_LP > tmin & time_LP < tmax);
        
        temp = data_x{ii,jj,1};
        fft_sig = (fft(temp));
        
        bins = find(freq_ax > 1e5);
        fft_sig(int32(bins)) = 0.0;
        fft_sig(npts - int32((bins))) = 0.0;
        inv_sig = ifft(fft_sig);

        filt_x = real(inv_sig);
        
        data_2x(ii,jj) = mean(filt_x(t_Langmuir));
        temp2 = data_y{ii,jj,1};
        data_2y(ii,jj) = mean(filt_y(t_Langmuir));
        temp3 = data_z{ii,jj,1};
        data_2z(ii,jj) = mean(filt_z(t_Langmuir));
        
        data_T1dx(ii,jj) = (filt_x(T1_ax));
    %     data_T1dy(1,ii) = data_y{1,ii,1}(T1_ax);
    %     data_T1dz(1,ii) = data_z{1,ii,1}(T1_ax);
        data_T2dx(ii,jj) = (filt_x(T2_ax));
    %     data_T2dy(1,ii) = data_y{1,ii,1}(T2_ax);
    %     data_T2dz(1,ii) = data_z{1,ii,1}(T2_ax);
        data_T3dx(ii,jj) = (filt_x(T3_ax));
    %     data_T3dy(1,ii) = data_y{1,ii,1}(T3_ax);
    %     data_T3dz(1,ii) = data_z{1,ii,1}(T3_ax);
        data_T4dx(ii,jj) = (filt_x(T4_ax));
    %     data_T4dy(1,ii) = data_y{1,ii,1}(T4_ax); 
    %     data_T4dz(1,ii) = data_z{1,ii,1}(T4_ax);
        
    end
end

%% Isat image plots
    
figure(4)
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

figure(5)
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














