% ----------------------------------------------------------------------- %
% Script to read and plot matlab (.mat) files for the coupled 1D 
% wave / transport solver. 
%
% rlbarnett c3149416 191106
% ----------------------------------------------------------------------- %

data = 1;
equilibrium = 0;
plots = 1;
para = 1;
perp = 0;

if data
    
    lapd_equib;
    clearvars -except n_init
    
    filepath = '/Volumes/DATA/LAPD/matlab/results_jsource_kyzero_eperpmix_mult07e5/';
    
    count = 1;

    for ii=106145

        filename = strcat(filepath, 'coupled_transport_', num2str(ii),'.mat');

        load(filename)
        
        if exist('rf_e','var')
            rf_ez = rf_e(1,3:3:3*npts);
            rf_ey = rf_e(1,2:3:3*npts);
            rf_ex = rf_e(1,1:3:3*npts);
        else
        end
        
        den_pert(count,:) = (n_init - n_new)./n_init;
        count = count + 1;
%         
%         x0 = 0;
%         y0 = 0;
%         width = 1300;
%         height = 600;
%         
% %         n = 8;
% %         CM = magma(n);
% 
%         col_purp = [0.28 0.06 0.47];
%         
%         cols = ["k","b","--b","r"];
%     
%         figure(1)
%         set(gcf,'Position',[x0 y0 width height],'color','w')
% %         subplot(2,1,1)
%         plot(nxax,den_pert,'k','linewidth',1.0,'DisplayName',...
%             ['{\itt} = ' num2str(round(double(ii)*dt/period)) ' T_{RF}'])
%         hold on
%         xlabel('Position (m)','interpreter','latex')
%         ylabel('$R_n (\times10^{-2}$)','interpreter','latex')
%         xlim([min(nxax) max(nxax)])
%         ylim([-0.02 0.02])
%         yticks([-0.02 -0.01 0.0 0.01 0.02])
%         yticklabels([-2.0 -1.0 0.0 1.0 2.0])
%         set(gca,'Fontsize',25)
%         legend('location','southwest')
%         legend('show')
% 
%         count = count+1;
        
    end
end

%%


x0 = 0;
y0 = 0;
width = 420*2;
height = 590*2;

figure(1)
set(gcf,'Units','points','Position',[x0 y0 width height],'Color','w')
subplot(5,1,1)
plot(nxax,den_pert(1,:),'k','linewidth',1.0)
hold on
ylabel('$R_n (\times10^{-2}$)','interpreter','latex')
xlim([min(nxax) max(nxax)])
xticks([])
ylim([-0.02 0.02])
yticks([-0.02 -0.01 0.0 0.01 0.02])
yticklabels([-2.0 -1.0 0.0 1.0 2.0])
set(gca,'Fontsize',20)
legend('$t = 10$ T$_{RF}$',...
    'interpreter','latex')

subplot(5,1,2)
plot(nxax,den_pert(2,:),'k','linewidth',1.0)
hold on
ylabel('$R_n (\times10^{-2}$)','interpreter','latex')
xlim([min(nxax) max(nxax)])
xticks([])
ylim([-0.02 0.02])
yticks([-0.02 -0.01 0.0 0.01 0.02])
yticklabels([-2.0 -1.0 0.0 1.0 2.0])
set(gca,'Fontsize',20)
legend('$t = 25$ T$_{RF}$',...
    'interpreter','latex')

subplot(5,1,3)
plot(nxax,den_pert(3,:),'k','linewidth',1.0)
hold on
ylabel('$R_n (\times10^{-2}$)','interpreter','latex')
xlim([min(nxax) max(nxax)])
xticks([])
ylim([-0.02 0.02])
yticks([-0.02 -0.01 0.0 0.01 0.02])
yticklabels([-2.0 -1.0 0.0 1.0 2.0])
set(gca,'Fontsize',20)
legend('$t = 50$ T$_{RF}$',...
    'interpreter','latex')

subplot(5,1,4)
plot(nxax,den_pert(4,:),'k','linewidth',1.0)
hold on
xticks([])
ylabel('$R_n (\times10^{-2}$)','interpreter','latex')
xlim([min(nxax) max(nxax)])
ylim([-0.02 0.02])
yticks([-0.02 -0.01 0.0 0.01 0.02])
yticklabels([-2.0 -1.0 0.0 1.0 2.0])
set(gca,'Fontsize',20)
legend('$t = 75$ T$_{RF}$',...
    'interpreter','latex')

subplot(5,1,5)
plot(nxax,den_pert(5,:),'k','linewidth',1.0)
hold on
xlabel('Position (m)','interpreter','latex')
ylabel('$R_n (\times10^{-2}$)','interpreter','latex')
xlim([min(nxax) max(nxax)])
ylim([-0.02 0.02])
yticks([-0.02 -0.01 0.0 0.01 0.02])
yticklabels([-2.0 -1.0 0.0 1.0 2.0])
set(gca,'Fontsize',20)
legend('$t = 100$ T$_{RF}$',...
    'interpreter','latex')

%%

% export_fig('/Volumes/DATA/thesis/figs/Rn_wtime_eperp_subplots.png',...
%     '-r300')