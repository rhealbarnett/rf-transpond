%-------------------------------------------------------------------------%
% read coupled model scan results
% rlbarnett c3149416 200630
%-------------------------------------------------------------------------%

filepath = '/Volumes/DATA/matlab/scan_results/powerscan_gaussian_eperp/';

%--
% string arrays for kx and density
kx_strarr = ["0+10i","0+15i","0+20i","0+25i"];
% kx_strarr = ["0+20i"];
% n_new_strarr = ["5e+16","7.925e+16","1.256e+17","1.9907e+17","3.155e+17","5e+17","8e+17","1.25e+18"];
n_new_strarr = ["5e+16","5e+17"];

%--
% cell array containing multipliers for the source term. 
%--
% kW_arr{1:5} is {100,200,400,800,1000}kW
%--
% kW_arr{}(1:4) is a string of the source multiplier for each kx
% e.g. kW_arr{1}(1:4) = ["mult for 10i","mult for 15i","mult for 20i"...]
% for 100kW. 
%--
kW_arr = cell(5,1);
kW_arr{1} = ["43100","55000","67200","81000"];
kW_arr{2} = ["60700","77600","95000","114000"];
kW_arr{3} = ["85800","109800","134300","161400"];
kW_arr{4} = ["122000","155700","190500","228600"];
kW_arr{5} = ["136200","174000","212500","255400"];
% kW_arr{1} = ["2363000"];
% kW_arr{2} = ["3346000"];
% kW_arr{3} = ["4730000"];
% kW_arr{4} = ["6683000"];
% kW_arr{5} = ["7473000"];

final_iter = '106145';

den_plot = cell(5,1);
n_arr = cell(5,1);

for ww=1:5
    for kk=1:4
        for jj=1
            if jj==7 || jj==8
                filename = strcat(filepath, 'coupled_',final_iter,'_',kW_arr{ww}(kk),...
                '_',kx_strarr(kk),'_',n_new_strarr(jj),'total.mat');
            else
%                 filename = strcat(filepath, 'coupled_',final_iter,'_',kW_arr{ww}(kk),...
%                 '_',kx_strarr(kk),'_',n_new_strarr(jj),'.mat');
                filename = strcat(filepath,'coupled_transport_106145.mat');
            end
            try
                load(filename);
                disp('Success!');
            catch LE
                fprintf('Nope: %s\n',LE.message);
            end
            fprintf('kx = %s, n_init = %s\n',[kx_strarr(kk),n_new_strarr(jj)])
%             fprintf('density perturbation %d\n',mean_pert)
%             den_plot{ww}(kk,jj) = mean_pert;
            n_arr{ww}(:,kk,jj) = n_new;
        end
    end
end


n_new_arr = [0.5e17,0.7925e17,1.256e17,1.9907e17,3.155e17,5.0e17,7.925e17,1.25e18];
kx_arr = [10i,15i,20i,25i];
power_arr = [100,200,400,800,1000];

%%


x0 = 0;
y0 = 0;
width = 1200;
height = 500;

figure(1)
set(gcf,'Position',[x0 y0 width height],'color','w')
plot(power_arr,cellfun(@(v)v(1),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(1),den_plot),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(3),den_plot),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(4),den_plot),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(5),den_plot),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(6),den_plot),'--+r','linewidth',1.5,'markersize',10)
% plot(power_arr,cellfun(@(v)v(7),den_plot),'-ok','linewidth',1.5,'markersize',10)
% plot(power_arr,cellfun(@(v)v(8),den_plot),'-ob','linewidth',1.5,'markersize',10)
% ylim([0.0 1.0e-4])
% yticks([0 0.5 1.0]*1e-4)
% yticklabels([{'0.0'} {'0.5'} {'1.0'}])
ylabel('$\overline{R_n}$ $\times10^{-4}$','interpreter','latex')
xlabel('Power (kW)','interpreter','latex')
legend('$n_{max} = 0.5e17$','$n_{max} = 0.8e17$','$n_{max} = 1.25e17$',...
    '$n_{max} = 2.0e17$','$n_{max} = 3.15e17$','$n_{max} = 5.0e17$',...
    '$n_{max} = 8.0e17$','$n_{max} = 1.25e18$','interpreter','latex',...
    'location','northwest','numcolumns',2)
% 
% export_fig('/Volumes/DATA/thesis/figs/nvp_deltaone_fixedscaling_epara.png',...
%     '-r300')

%%

x0 = 0;
y0 = 0;
width = 420*2;
height = 590*2;

figure(2)
set(gcf,'Units','points','Position',[x0 y0 width height],'Color','w')
subplot(4,1,1)
plot(power_arr,cellfun(@(v)v(1),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(5),den_plot),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(9),den_plot),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(13),den_plot),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(17),den_plot),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(21),den_plot),'--+r','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(25),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(29),den_plot),'-ob','linewidth',1.5,'markersize',10)
xticks([])
ylim([0.0 0.08])
yticks([0.0 0.04 0.08])
yticklabels([{'0.0'} {'4.0'} {'8.0'}])
ylabel('$\overline{R_n}$ $\times10^{-2}$','interpreter','latex')
text(0.005,0.98,'(a) $k_x = 10i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,2)
plot(power_arr,cellfun(@(v)v(2),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(6),den_plot),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(10),den_plot),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(14),den_plot),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(18),den_plot),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(22),den_plot),'--+r','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(26),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(30),den_plot),'-ob','linewidth',1.5,'markersize',10)
xticks([])
ylim([0.0 0.08])
yticks([0.0 0.04 0.08])
yticklabels([{'0.0'} {'4.0'} {'8.0'}])
ylabel('$\overline{R_n}$ $\times10^{-2}$','interpreter','latex')
text(0.005,0.98,'(b) $k_x = 15i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,3)
plot(power_arr,cellfun(@(v)v(3),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(7),den_plot),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(11),den_plot),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(15),den_plot),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(19),den_plot),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(23),den_plot),'--+r','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(27),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(31),den_plot),'-ob','linewidth',1.5,'markersize',10)
xticks([])
ylim([0.0 0.08])
yticks([0.0 0.04 0.08])
yticklabels([{'0.0'} {'4.0'} {'8.0'}])
ylabel('$\overline{R_n}$ $\times10^{-2}$','interpreter','latex')
text(0.005,0.98,'(c) $k_x = 20i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,4)
plot(power_arr,cellfun(@(v)v(4),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(8),den_plot),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(12),den_plot),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(16),den_plot),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(20),den_plot),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(24),den_plot),'--+r','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(28),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(32),den_plot),'-ob','linewidth',1.5,'markersize',10)
xlabel('Power (kW)','interpreter','latex')
ylim([0.0 0.08])
yticks([0.0 0.04 0.08])
yticklabels([{'0.0'} {'4.0'} {'8.0'}])
ylabel('$\overline{R_n}$ $\times10^{-2}$','interpreter','latex')
text(0.005,0.98,'(d) $k_x = 25i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
legend({'$n_{max} = 0.5e17$','$n_{max} = 0.8e17$','$n_{max} = 1.25e17$'...
    '$n_{max} = 2.0e17$','$n_{max} = 3.15e17$','$n_{max} = 5.0e17$',...
    '$n_{max} = 8.0e17$','$n_{max} = 1.25e18$'},'interpreter','latex',...
    'location','northeast','NumColumns',2)
hold off

% export_fig('/Volumes/DATA/thesis/figs/Rn_vs_power_eperp_subplots_extran.png',...
%     '-r300')

%%

den_pert_ratio = cell(5,1);
den_pert_diff = cell(5,1);
den_pert_arr = cell(5,1);
n_init_arr = cell(5,1);
tfmin_arr = cell(5,1);
tfmax_arr = cell(5,1);
den_pert_min = cell(5,1);
den_pert_max = cell(5,1);

for ww=1:5
    for kk=1:4
        for jj=1:2
            ww
            kk
            jj
            if jj==7
                filename = strcat(filepath, 'coupled_',final_iter,'_',kW_arr{ww}(kk),...
                '_',kx_strarr(kk),'_',n_new_strarr(jj),'total.mat');
            else 
                filename = strcat(filepath, 'coupled_',final_iter,'_',kW_arr{ww}(kk),...
                '_',kx_strarr(kk),'_',n_new_strarr(jj),'.mat');
            end
            Nmax = str2double(n_new_strarr(jj));
%             kx = str2double(kx_strarr(kk));
%             lapd_equib;
%             clearvars -except n_init filename filepath;
            equib = load('equil_transport_input.mat','n_new');
            fact = Nmax/max(equib.n_new);
            equib.n_new = equib.n_new*fact;
            n_init = equib.n_new;
            clear equib.n_new
            try
                load(filename);
                disp('Success!');
            catch LE
                fprintf('Nope: %s\n',LE.message);
            end
%             clear mean_pert
            rf_source_varax = interp1(zax,rf_source,nxax,'linear');
            rf_source_varax(isnan(rf_source_varax)) = 0;
            density_perturbation = (n_init - n_new) ./ n_init;
            ax = find(nxax >= -0.5 & nxax <= 0.5);
%             wind_np = length(ax);
%             tfmin = islocalmin(density_perturbation(ax));
            tfmax = islocalmax(density_perturbation(ax));
            tfmin = islocalmin(density_perturbation(ax));
%             tfmin(105:106) = 0;
%             tfmax(105:106) = 0;
            for gg=floor(length(ax)/2):length(ax)
%                 if tfmin(gg)
%                     tfmin(gg+1:end) = 0;
%                     tfmin(1:length(ax)-gg) = 0;
%                 end
                if tfmax(gg)
                    tfmax(gg+1:end) = 0;
                    tfmax(1:length(ax)-gg) = 0;
                end
            end
%             max_ind = ax(tfmax);
%             tfmin = islocalmin(density_perturbation(max_ind(1):max_ind(2)));
%             min_ind = find(tfmin==1);
            for gg=1:(length(ax))
                if sum(tfmax(1:gg))==1
                    tfmin(gg) = 0;
%                     tfmin(length(ax)-gg+1) = 0;
%                 elseif tfmin(gg) && sum(tfmax(1:gg))==2
%                     tfmin(gg+1:end) = 0;
%                     tfmin(1:length(ax)-gg) = 0;
                end
            end
            for gg=floor(length(ax)/2):length(ax)
                if sum(tfmin(floor(length(ax)/2):gg))==1
                    tfmin(gg+1:end) = 0;
                    tfmin(tfmin(1:length(ax)-gg)) = 0;
                end
                
            end
%             plot(nxax(ax),density_perturbation(ax),'k')
%             hold on
%             plot(nxax(ax(tfmin)),density_perturbation(ax(tfmin)),'*b')
%             plot(nxax(ax(tfmax)),density_perturbation(ax(tfmax)),'*r')
%             close 1
            min_mean = abs(mean(density_perturbation(ax(tfmin))));
            max_mean = abs(mean(density_perturbation(ax(tfmax))));
            tfmin_arr{ww}(kk,jj,:) = tfmin;
            tfmax_arr{ww}(kk,jj,:) = tfmax;
            den_pert_ratio{ww}(kk,jj) = min_mean/max_mean;
            den_pert_diff{ww}(kk,jj) = max_mean - min_mean;
            den_pert_min{ww}(kk,jj) = mean(density_perturbation(ax(tfmin)));
            den_pert_max{ww}(kk,jj) = mean(density_perturbation(ax(tfmax)));
            mean(density_perturbation(ax(tfmin)))
            mean(density_perturbation(ax(tfmax)))
%             indmin = find(tfmin==1);
%             indmax = find(tfmax==1);
%             uni_ax = linspace(min(nxax(ax)),max(nxax(ax)),wind_np);
%             wind = hann(wind_np);
%             pert_wind = interp1(uni_ax,wind,nxax(ax),'linear');
%             density_wind = pert_wind.*density_perturbation(ax);
%             mean_pert = mean(density_wind);
%             fprintf('kx = %s, n_init = %s\n',[kx_strarr(kk),n_new_strarr(jj)])
%             fprintf('density perturbation %d\n',mean_pert)
            den_plot{ww}(kk,jj) = mean_pert;
            den_pert_arr{ww}(:,kk,jj) = density_perturbation;
            n_init_arr{ww}(:,kk,jj) = n_init;
            n_arr{ww}(:,kk,jj) = n_new;
        end
    end
end

%%

n_new_arr = [5.0e16,5.0e17];
kx_arr = [10i,15i,20i,25i];
power_arr = [100,200,400,800,1000];

x0 = 0;
y0 = 0;
width = 1200;
height = 500;

figure(3)
set(gcf,'Position',[x0 y0 width height],'color','w')
plot(power_arr,cellfun(@(v)v(1),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(5),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(1),den_pert_min),'-xb','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(1),den_pert_max),'-xr','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(5),den_pert_min),'-ob','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(5),den_pert_max),'-or','linewidth',1.,'markersize',10)
% plot(power_arr,cellfun(@(v)v(5),den_pert_arr),'-xb','linewidth',1.5,'markersize',10)
% plot(power_arr,cellfun(@(v)v(11),den_pert_arr),'-xr','linewidth',1.5,'markersize',10)
% plot(power_arr,cellfun(@(v)v(15),den_pert_arr),'--+k','linewidth',1.5,'markersize',10)
% plot(power_arr,cellfun(@(v)v(19),den_pert_arr),'--+b','linewidth',1.5,'markersize',10)
% plot(power_arr,cellfun(@(v)v(23),den_pert_arr),'--+r','linewidth',1.5,'markersize',10)
% ylim([-0.04 0.1])
ylabel('$\overline{R_n}$','interpreter','latex')
xlabel('Power (kW)','interpreter','latex')
% legend('$n_{max} = 0.5e17$','$n_{max} = 0.8e17$','$n_{max} = 1.25e17$',...
%     '$n_{max} = 2.0e17$','$n_{max} = 3.15e17$','$n_{max} = 5.0e17$','interpreter','latex',...
%     'location','northwest')
% legend({'$n_{max} = 1.25e17$'...
%     '$n_{max} = 1.25e18$'},'interpreter','latex',...
%     'location','northeast')

% export_fig('/Volumes/DATA/thesis/figs/nvp_kx25i_diff_fixedyscale.png',...
%     '-r300')

%%

x0 = 0;
y0 = 0;
width = 420*2;
height = 600*2;

figure(4)
set(gcf,'Units','points','Position',[x0 y0 width height],'Color','w')
subplot(4,1,1)
plot(power_arr,cellfun(@(v)v(1),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(5),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(1),den_pert_min),'-xb','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(1),den_pert_max),'-xr','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(5),den_pert_min),'-ob','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(5),den_pert_max),'-or','linewidth',1.,'markersize',10)
xticks([])
ylim([-0.23 0.23])
yticks([-0.2 0.0 0.2])
yticklabels([{'-2.0'} {'1.0'} {'2.0'}])
ylabel('$\overline{R_n}$ $\times10^{-1}$','interpreter','latex')
text(0.005,0.98,'(a) $k_x = 10i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,2)
plot(power_arr,cellfun(@(v)v(2),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(6),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(2),den_pert_min),'-xb','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(2),den_pert_max),'-xr','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(6),den_pert_min),'-ob','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(6),den_pert_max),'-or','linewidth',1.,'markersize',10)
xticks([])
ylim([-0.23 0.23])
yticks([-0.2 0.0 0.2])
yticklabels([{'-2.0'} {'1.0'} {'2.0'}])
ylabel('$\overline{R_n}$ $\times10^{-1}$','interpreter','latex')
text(0.005,0.98,'(b) $k_x = 15i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,3)
plot(power_arr,cellfun(@(v)v(3),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(7),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(3),den_pert_min),'-xb','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(3),den_pert_max),'-xr','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(7),den_pert_min),'-ob','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(7),den_pert_max),'-or','linewidth',1.,'markersize',10)
xticks([])
ylim([-0.23 0.23])
yticks([-0.2 0.0 0.2])
yticklabels([{'-2.0'} {'1.0'} {'2.0'}])
ylabel('$\overline{R_n}$ $\times10^{-1}$','interpreter','latex')
text(0.005,0.98,'(c) $k_x = 20i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,4)
h(1) = plot(power_arr,cellfun(@(v)v(4),den_plot),'-xk','linewidth',1.5,'markersize',10);
hold on
h(2) = plot(power_arr,cellfun(@(v)v(8),den_plot),'-ok','linewidth',1.5,'markersize',10);
h(3) = plot(power_arr,cellfun(@(v)v(4),den_pert_min),'-xb','linewidth',1.,'markersize',10);
h(4) = plot(power_arr,cellfun(@(v)v(4),den_pert_max),'-xr','linewidth',1.,'markersize',10);
h(5) = plot(power_arr,cellfun(@(v)v(8),den_pert_min),'-ob','linewidth',1.,'markersize',10);
h(6) = plot(power_arr,cellfun(@(v)v(8),den_pert_max),'-or','linewidth',1.,'markersize',10);
xlabel('Power (kW)','interpreter','latex')
ylim([-0.23 0.23])
yticks([-0.2 0.0 0.2])
yticklabels([{'-2.0'} {'1.0'} {'2.0'}])
ylabel('$\overline{R_n}$ $\times10^{-1}$','interpreter','latex')
text(0.005,0.98,'(d) $k_x = 25i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
legend(h([1,2]),{'$n_{max} = 1.25e17$'...
    '$n_{max} = 1.25e18$'},'interpreter','latex',...
    'location','northeast','NumColumns',2)
hold off

% export_fig('/Volumes/DATA/thesis/figs/meanRn_pm_extrema.png',...
%     '-r300')

%%

x0 = 0;
y0 = 0;
width = 1200;
height = 500;

figure(4)
set(gcf,'Position',[x0 y0 width height],'color','w')
plot(nxax,den_pert_arr{1}(:,3),'k','linewidth',0.9,'markersize',10)
hold on
plot(nxax,den_pert_arr{2}(:,3),'b','linewidth',0.9,'markersize',10)
plot(nxax,den_pert_arr{3}(:,3),'r','linewidth',0.9,'markersize',10)
plot(nxax,den_pert_arr{4}(:,3),'--k','linewidth',0.9,'markersize',10)
plot(nxax,den_pert_arr{5}(:,3),'--r','linewidth',0.9,'markersize',10)
xlim([-.5 .5])
ylim([-6.0e-2 6.0e-2])
yticks([-6.0e-2 -3.0e-2 0.0 3.0e-2 6.0e-2])
yticklabels([-6.0 -3.0 0.0 3.0 6.0])
ylabel('$R_n$ $(\times10^{-2})$','interpreter','latex')
xlabel('Position (m)','interpreter','latex')
legend('$P = 100$ kW','$P = 200$ kW','$P = 400$ kW',...
    '$P = 800$ kW','$P = 1000$ kW','interpreter','latex',...
    'location','south')

% export_fig('/Volumes/DATA/thesis/figs/Rnvsn_powerscan_eperp.png',...
%     '-r300')

%%

x0 = 0;
y0 = 0;
width = 420*2;
height = 590*2;

figure(4)
set(gcf,'Units','points','Position',[x0 y0 width height],'Color','w')
subplot(4,1,1)
plot(power_arr,cellfun(@(v)v(1),den_pert_diff),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(5),den_pert_diff),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(9),den_pert_diff),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(13),den_pert_diff),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(17),den_pert_diff),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(21),den_pert_diff),'--+r','linewidth',1.5,'markersize',10)
xticks([])
ylim([-2.5e-5 1.6e-4])
yticks([0 0.5 1 1.5]*1e-4)
yticklabels([{'0.0'} {'0.5'} {'1.0'} {'1.5'}])
ylabel('$\overline{R_n}$ $\times10^{-4}$','interpreter','latex')
text(0.005,0.98,'(a) $k_x = 10i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,2)
plot(power_arr,cellfun(@(v)v(2),den_pert_diff),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(6),den_pert_diff),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(10),den_pert_diff),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(14),den_pert_diff),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(18),den_pert_diff),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(22),den_pert_diff),'--+r','linewidth',1.5,'markersize',10)
xticks([])
ylim([-2.5e-5 1.6e-4])
yticks([0 0.5 1 1.5]*1e-4)
yticklabels([{'0.0'} {'0.5'} {'1.0'} {'1.5'}])
ylabel('$\overline{R_n}$ $\times10^{-4}$','interpreter','latex')
text(0.005,0.98,'(b) $k_x = 15i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,3)
plot(power_arr,cellfun(@(v)v(3),den_pert_diff),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(7),den_pert_diff),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(11),den_pert_diff),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(15),den_pert_diff),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(19),den_pert_diff),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(23),den_pert_diff),'--+r','linewidth',1.5,'markersize',10)
xticks([])
ylim([-2.5e-5 1.6e-4])
yticks([0 0.5 1 1.5]*1e-4)
yticklabels([{'0.0'} {'0.5'} {'1.0'} {'1.5'}])
ylabel('$\overline{R_n}$ $\times10^{-4}$','interpreter','latex')
text(0.005,0.98,'(c) $k_x = 20i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
hold off

subplot(4,1,4)
plot(power_arr,cellfun(@(v)v(4),den_pert_diff),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(8),den_pert_diff),'-xb','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(12),den_pert_diff),'-xr','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(16),den_pert_diff),'--+k','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(20),den_pert_diff),'--+b','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(24),den_pert_diff),'--+r','linewidth',1.5,'markersize',10)
xlabel('Power (kW)','interpreter','latex')
ylim([-2.5e-5 1.6e-4])
yticks([0 0.5 1 1.5]*1e-4)
yticklabels([{'0.0'} {'0.5'} {'1.0'} {'1.5'}])
ylabel('$\overline{R_n}$ $\times10^{-4}$','interpreter','latex')
text(0.005,0.98,'(d) $k_x = 25i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')
legend({'$n_{max} = 0.5e17$','$n_{max} = 0.8e17$','$n_{max} = 1.25e17$'...
    '$n_{max} = 2.0e17$','$n_{max} = 3.15e17$','$n_{max} = 5.0e17$'},'interpreter','latex',...
    'location','northeast','NumColumns',2)
hold off

% export_fig('/Volumes/DATA/thesis/figs/nvp_diff_subplots_epara.png',...
%     '-r300')

%%

figure(5)
set(gcf,'Position',[x0 y0 width height],'color','w')
subplot(4,1,1)
plot(nxax(ax),den_pert_arr{4}(ax,1,2),'k','linewidth',1.)
hold on
plot(nxax(ax(tfmin_arr{4}(1,2,:))),den_pert_arr{4}(ax(tfmin_arr{4}(1,2,:)),1,2),'xb','markersize',10,'linewidth',2.)
plot(nxax(ax(tfmax_arr{4}(1,2,:))),den_pert_arr{4}(ax(tfmax_arr{4}(1,2,:)),1,2),'xr','markersize',10,'linewidth',2.)
xticks([])
ylim([-0.1 0.2])
yticks([-0.1 0.0 0.1 0.2])
yticklabels([{'-1.0'} {'0.0'} {'1.0'} {'2.0'}])
% xlabel('Position (m)','interpreter','latex')
ylabel('$R_n$ $\times10^{-1}$','interpreter','latex')
text(0.005,0.98,'(a) $k_x = 10i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')

subplot(4,1,2)
plot(nxax(ax),den_pert_arr{4}(ax,2,2),'k','linewidth',1.)
hold on
plot(nxax(ax(tfmin_arr{4}(2,2,:))),den_pert_arr{4}(ax(tfmin_arr{4}(2,2,:)),2,2),'xb','markersize',10,'linewidth',2.)
plot(nxax(ax(tfmax_arr{4}(2,2,:))),den_pert_arr{4}(ax(tfmax_arr{4}(2,2,:)),2,2),'xr','markersize',10,'linewidth',2.)
xticks([])
ylim([-0.1 0.2])
yticks([-0.1 0.0 0.1 0.2])
yticklabels([{'-1.0'} {'0.0'} {'1.0'} {'2.0'}])
% xlabel('Position (m)','interpreter','latex')
ylabel('$R_n$ $\times10^{-1}$','interpreter','latex')
text(0.005,0.98,'(b) $k_x = 15i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')

subplot(4,1,3)
plot(nxax(ax),den_pert_arr{4}(ax,3,2),'k','linewidth',1.)
hold on
plot(nxax(ax(tfmin_arr{4}(3,2,:))),den_pert_arr{4}(ax(tfmin_arr{4}(3,2,:)),3,2),'xb','markersize',10,'linewidth',2.)
plot(nxax(ax(tfmax_arr{4}(3,2,:))),den_pert_arr{4}(ax(tfmax_arr{4}(3,2,:)),3,2),'xr','markersize',10,'linewidth',2.)
xticks([])
ylim([-0.1 0.2])
yticks([-0.1 0.0 0.1 0.2])
yticklabels([{'-1.0'} {'0.0'} {'1.0'} {'2.0'}])
% xlabel('Position (m)','interpreter','latex')
ylabel('$R_n$ $\times10^{-1}$','interpreter','latex')
text(0.005,0.98,'(c) $k_x = 20i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')

subplot(4,1,4)
plot(nxax(ax),den_pert_arr{4}(ax,4,2),'k','linewidth',1.)
hold on
plot(nxax(ax(tfmin_arr{4}(4,2,:))),den_pert_arr{4}(ax(tfmin_arr{4}(4,2,:)),4,2),'xb','markersize',10,'linewidth',2.)
plot(nxax(ax(tfmax_arr{4}(4,2,:))),den_pert_arr{4}(ax(tfmax_arr{4}(4,2,:)),4,2),'xr','markersize',10,'linewidth',2.)
ylim([-0.1 0.2])
yticks([-0.1 0.0 0.1 0.2])
yticklabels([{'-1.0'} {'0.0'} {'1.0'} {'2.0'}])
xlabel('Position (m)','interpreter','latex')
ylabel('$R_n$ $\times10^{-1}$','interpreter','latex')
text(0.005,0.98,'(d) $k_x = 25i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',20,...
                'color','black','interpreter','latex')

% export_fig('/Volumes/DATA/thesis/figs/npert_extrema_n1e18_800kW.png',...
%     '-r300')

%%

x0 = 0;
y0 = 0;
width = 1200;
height = 500;

figure(6)
set(gcf,'Position',[x0 y0 width height],'color','w')
plot(nxax(ax),den_pert_arr{1}(ax,1,1),'k')
hold on
plot(nxax(ax),den_pert_arr{2}(ax,1,1),'b')
plot(nxax(ax),den_pert_arr{3}(ax,1,1),'r')
plot(nxax(ax),den_pert_arr{4}(ax,1,1),'--k')
plot(nxax(ax),den_pert_arr{5}(ax,1,1),'--b')
legend('100 kW','200 kW','400 kW','800 kW','1 MW','interpreter','latex')
plot(nxax(ax(tfmin_arr{1}(1,1,:))),den_pert_arr{1}(ax(tfmin_arr{1}(1,1,:)),1,1),'*k','HandleVisibility','off')
plot(nxax(ax(tfmax_arr{1}(1,1,:))),den_pert_arr{1}(ax(tfmax_arr{1}(1,1,:)),1,1),'*k','HandleVisibility','off')
plot(nxax(ax(tfmin_arr{2}(1,1,:))),den_pert_arr{2}(ax(tfmin_arr{2}(1,1,:)),1,1),'*b','HandleVisibility','off')
plot(nxax(ax(tfmax_arr{2}(1,1,:))),den_pert_arr{2}(ax(tfmax_arr{2}(1,1,:)),1,1),'*b','HandleVisibility','off')
plot(nxax(ax(tfmin_arr{3}(1,1,:))),den_pert_arr{3}(ax(tfmin_arr{3}(1,1,:)),1,1),'*r','HandleVisibility','off')
plot(nxax(ax(tfmax_arr{3}(1,1,:))),den_pert_arr{3}(ax(tfmax_arr{3}(1,1,:)),1,1),'*r','HandleVisibility','off')
plot(nxax(ax(tfmin_arr{4}(1,1,:))),den_pert_arr{4}(ax(tfmin_arr{4}(1,1,:)),1,1),'ok','HandleVisibility','off')
plot(nxax(ax(tfmax_arr{4}(1,1,:))),den_pert_arr{4}(ax(tfmax_arr{4}(1,1,:)),1,1),'ok','HandleVisibility','off')
plot(nxax(ax(tfmin_arr{5}(1,1,:))),den_pert_arr{5}(ax(tfmin_arr{5}(1,1,:)),1,1),'ob','HandleVisibility','off')
plot(nxax(ax(tfmax_arr{5}(1,1,:))),den_pert_arr{5}(ax(tfmax_arr{5}(1,1,:)),1,1),'ob','HandleVisibility','off')
xlabel('Position (m)','interpreter','latex')
ylabel('$R_n$','interpreter','latex')

% export_fig('/Volumes/DATA/thesis/figs/Rn_extremalocs_perp_kx10_n5e16.png',...
%     '-r300')

%%

x0 = 0;
y0 = 0;
width = 1500;
height = 700;

figure(4)
set(gcf,'Units','points','Position',[x0 y0 width height],'Color','w')
subplot(2,2,1)
plot(power_arr,cellfun(@(v)v(1),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(5),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(1),den_pert_min),'-xb','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(1),den_pert_max),'-xr','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(5),den_pert_min),'-ob','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(5),den_pert_max),'-or','linewidth',1.,'markersize',10)
xticks([])
ylim([-0.3 0.3])
xlim([100,1000])
yticks([-0.3 0.0 0.3])
yticklabels([{'-3.0'} {'0.0'} {'3.0'}])
ylabel('$\overline{R_n}$ $\times10^{-1}$','interpreter','latex')
set(gca,'Fontsize',30,'FontName','CMU Serif')
text(0.005,0.98,'$k_x = 10i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','black','interpreter','latex')
hold off

subplot(2,2,3)
plot(power_arr,cellfun(@(v)v(2),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(6),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(2),den_pert_min),'-xb','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(2),den_pert_max),'-xr','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(6),den_pert_min),'-ob','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(6),den_pert_max),'-or','linewidth',1.,'markersize',10)
% xticks([])
xlabel('Power (kW)','interpreter','latex')
ylim([-0.3 0.3])
xlim([100,1000])
yticks([-0.3 0.0 0.3])
yticklabels([{'-3.0'} {'0.0'} {'3.0'}])
ylabel('$\overline{R_n}$ $\times10^{-1}$','interpreter','latex')
set(gca,'Fontsize',30,'FontName','CMU Serif')
text(0.005,0.98,'$k_x = 15i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','black','interpreter','latex')
hold off

subplot(2,2,2)
plot(power_arr,cellfun(@(v)v(3),den_plot),'-xk','linewidth',1.5,'markersize',10)
hold on
plot(power_arr,cellfun(@(v)v(7),den_plot),'-ok','linewidth',1.5,'markersize',10)
plot(power_arr,cellfun(@(v)v(3),den_pert_min),'-xb','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(3),den_pert_max),'-xr','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(7),den_pert_min),'-ob','linewidth',1.,'markersize',10)
plot(power_arr,cellfun(@(v)v(7),den_pert_max),'-or','linewidth',1.,'markersize',10)
xticks([])
ylim([-0.3 0.3])
xlim([100,1000])
% yticks([-0.2 0.0 0.2])
% yticklabels([{'-2.0'} {'1.0'} {'2.0'}])
% ylabel('$\overline{R_n}$ $\times10^{-1}$','interpreter','latex')
set(gca,'Fontsize',30,'FontName','CMU Serif')
yticks([])
text(0.005,0.98,'$k_x = 20i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','black','interpreter','latex')
hold off

subplot(2,2,4)
h(1) = plot(power_arr,cellfun(@(v)v(4),den_plot),'-xk','linewidth',1.5,'markersize',10);
hold on
h(2) = plot(power_arr,cellfun(@(v)v(8),den_plot),'-ok','linewidth',1.5,'markersize',10);
h(3) = plot(power_arr,cellfun(@(v)v(4),den_pert_min),'-xb','linewidth',1.,'markersize',10);
h(4) = plot(power_arr,cellfun(@(v)v(4),den_pert_max),'-xr','linewidth',1.,'markersize',10);
h(5) = plot(power_arr,cellfun(@(v)v(8),den_pert_min),'-ob','linewidth',1.,'markersize',10);
h(6) = plot(power_arr,cellfun(@(v)v(8),den_pert_max),'-or','linewidth',1.,'markersize',10);
xlabel('Power (kW)','interpreter','latex')
ylim([-0.3 0.3])
xlim([100,1000])
yticks([])
% yticks([-0.2 0.0 0.2])
% yticklabels([{'-2.0'} {'1.0'} {'2.0'}])
% ylabel('$\overline{R_n}$ $\times10^{-1}$','interpreter','latex')
set(gca,'Fontsize',30,'FontName','CMU Serif')
text(0.005,0.98,'$k_x = 25i$ m$^{-1}$','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
                'color','black','interpreter','latex')
legend(h([1,2]),{'$n_{max} = 5.0e16$'...
    '$n_{max} = 5.0e17$'},'interpreter','latex',...
    'location','northeast','NumColumns',2)
hold off
