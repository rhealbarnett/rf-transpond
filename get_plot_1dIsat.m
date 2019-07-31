%% Call in xz plane Isat data file

filename = '/Users/rhealbarnett/Documents/LAPD/RF_April22/11_Mach_p30xyz_12kV.hdf5';


%% get Mach probe data: check the LAPD_201904_documentation ---------------- %%

probename = 'Mach';
channels = ["Mach1", "Mach2", "Mach3", "Mach4","Mach5","Mach6","Antenna current"];


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

x = linspace(minx,maxx,Nx);
y = linspace(minz,maxz,Nz);

% attenuation (V) and resistance (omhs), needed to calibrate the signal.
atten = 3.16;
res = 25;

% time values taken 'safely' within the rf application.
rfmin = 0.7;
rfmax = 1.2;

% time values before the rf appliation.
tmin = 0;
tmax = 0.4;


%% get data ---------------------------------------------------------------- %

[actx_LP, acty_LP, time_LP, data_LP] = get_probe(filename, probename, channels, x, y, 1);

callib = atten/res;

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

z = y;

%% Calculate average of each probe for density proxy from Isat


for ii= 1:numel(x)
    for jj= 1:numel(y)
        data_x{ii,jj,1} = (data_LP{ii,jj,1} + data_LP{ii,jj,2})./2.0;
        data_y{ii,jj,1} = (data_LP{ii,jj,3} + data_LP{ii,jj,4})./2.0;
        data_z{ii,jj,1} = (data_LP{ii,jj,5} + data_LP{ii,jj,6})./2.0;
    end
end

%% Filter the data to remove the RF signal. 

for ii=1:numel(y)
    
    data_x{1,ii,1} = sgolayfilt(data_x{1,ii,1},3,1201);
    data_y{1,ii,1} = sgolayfilt(data_y{1,ii,1},3,1201);
    data_z{1,ii,1} = sgolayfilt(data_z{1,ii,1},3,1201);
    
end


%% For comparison to the 1D, parallel transport model:
% consider a line in the x-z plane, along z and as close to the antenna as
% possible (x=-11cm). 
% Average the Isat measurement before RF, ensuring to stop well before the
% RF starts. After the RF starts, excluding the transient period, look at
% snapshots at approximately 25, 50, 75 and 100 periods. Times are in ms.
% Approximate RF start is at 0.5008ms, end of transient behaviour approx
% 0.5348ms. 

f = 2.38e6;
T = 1.0/f;
T_ms = T*1.0e3;

TRF = 0.538;
T1 = 25*T_ms;
T2 = 50*T_ms;
T3 = 75*T_ms;
T4 = 100*T_ms;

T1_ax = find(time_LP<=(TRF+T1),1,'last');
T2_ax = find(time_LP<=(TRF+T2),1,'last');
T3_ax = find(time_LP<=(TRF+T3),1,'last');
T4_ax = find(time_LP<=(TRF+T4),1,'last');

for ii=1:length(z)
    data_T1dx(1,ii) = data_x{1,ii,1}(T1_ax);
    data_T1dy(1,ii) = data_y{1,ii,1}(T1_ax);
    data_T1dz(1,ii) = data_z{1,ii,1}(T1_ax);
    data_T2dx(1,ii) = data_x{1,ii,1}(T2_ax);
    data_T2dy(1,ii) = data_y{1,ii,1}(T2_ax);
    data_T2dz(1,ii) = data_z{1,ii,1}(T2_ax);
    data_T3dx(1,ii) = data_x{1,ii,1}(T3_ax);
    data_T3dy(1,ii) = data_y{1,ii,1}(T3_ax);
    data_T3dz(1,ii) = data_z{1,ii,1}(T3_ax);
    data_T4dx(1,ii) = data_x{1,ii,1}(T4_ax);
    data_T4dy(1,ii) = data_y{1,ii,1}(T4_ax);
    data_T4dz(1,ii) = data_z{1,ii,1}(T4_ax);
end

%%
% Use mean of the data values over the time period before RF is turned on.

for ii= 1:numel(x)
    for jj= 1:numel(z)
        
        t_Langmuir = find(time_LP > tmin & time_LP < tmax);
        
        temp1 = data_x{ii,jj,1};
        mach11bf(ii,jj) = mean(temp1(t_Langmuir));
        temp2 = data_y{ii,jj,1};
        mach12bf(ii,jj) = mean(temp2(t_Langmuir));
        temp3 = data_z{ii,jj,1};
        mach13bf(ii,jj) = mean(temp3(t_Langmuir));

        
    end
end

%%

figure(1)
plot(z,mach11bf(1,:),'.-')

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

plot(3.0*ones(1,Nz),linspace(0.006,0.02,length(data_T1dx)),'--k')
plot(-3.0*ones(1,Nz),linspace(0.006,0.02,length(data_T1dx)),'--k')

plot(5.0*ones(1,Nz),linspace(0.006,0.02,length(data_T1dx)),'--r')
plot(-5.0*ones(1,Nz),linspace(0.006,0.02,length(data_T1dx)),'--r')

hold off














