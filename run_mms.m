%-------------------------------------------------------------------------%
% run_mms.m function .. rhealbarnett 201028
%
% Runs MMS tests for the coupled transport equations. 
% 
% INPUTS:
% SS, TD : 0 or 1 to run steady state or time dependent test.
% mms_plots : 0 or 1 for plots. 
% save_file : 0 or 1 to not save/save MMS output. 
%
% OUTPUTS:
% linf(two)_arr : 5 X 2 matrix of L infinity (L^2) errors for each dx 
%            (row) and n (linf(two)_arr(:,1)), v (linf(two)_arr(:,2)). 
% ratio_inf(two)_arr : 4 X 2 matrix of error ratio for each dx (row) and n 
%            ratio_inf(two)_arr(:,1)), v (ratio_inf(two)_arr(:,2)). 
%            (see eqn 24 in Salari & Knupp, 2000)
% oo_inf(two)_arr : 4 X 2 matrix of order of accuracy for each dx (row) 
%            and n (oo_inf(two)_arr(:,1)), v (oo_inf(two)_arr(:,2)). 
%            (see eqn 25 in Salari & Knupp, 2000)
%            
%-------------------------------------------------------------------------%

function [linf_arr,ltwo_arr,ratio_inf_arr,ratio_two_arr,...
    oo_inf_arr,oo_two_arr] = run_mms(SS,TD,mms_plots,save_file)

    npts = 512;

    if SS

        dt = 1.0;
        tol = 1.0e-14;

    elseif TD

        dt = 8.0e-6;
        tmin = 0.0;
        tmax = 8.0e-4;

    end

    transport_mms;

    for kk=1:5

        if SS

            nmax = 1.0e5;
            save_iter = round(nmax/5);

            fprintf('iteration %d\n', kk)
            fprintf('Number of points %d\n', npts)

            fprintf('Minimum dx %d\n', min(ndx))
            fprintf('Maximum dx %d\n', max(ndx))

            transport_1d

            npts_arr(1,kk) = npts;
            npts = npts/2;

        elseif TD

            nmax = round(tmax/dt);
            save_iter = round(nmax/5);

            fprintf('iteration %d\n', kk)
            fprintf('Time step %d\n', dt)
            fprintf('nmax, tmax %d, %d\n', [nmax tmax])

            transport_1d

            dt_arr(1,kk) = dt;
            dt = dt/2.0;

        end

        ltwon_arr(1,kk) = l_twon;
        linfn_arr(1,kk) = l_infn;
        ltwou_arr(1,kk) = l_twou;
        linfu_arr(1,kk) = l_infu;

        if save_file

            if SS
                str = 'SS';
                mms_errors.npts_arr = npts_arr;
            elseif TD
                str = 'TD';
                mms_errors.dt_arr = dt_arr;
            end

            [status,git_hash] = system('git rev-parse HEAD');
            s1 = '# Created from matlab git hash ';
            s2 = git_hash;
            header = [s1 s2];

            mms_errors.header = header;
            mms_errors.ltwon_arr = ltwon_arr;
            mms_errors.linfn_arr = linfn_arr;
            mms_errors.ltwou_arr = ltwou_arr;
            mms_errors.linfu_arr = linfu_arr;

        end

    end


    ratio_infu = linfu_arr(1:kk-1)./linfu_arr(2:kk);
    ratio_twou = ltwou_arr(1:kk-1)./ltwou_arr(2:kk);
    ratio_infn = linfn_arr(1:kk-1)./linfn_arr(2:kk);
    ratio_twon = ltwon_arr(1:kk-1)./ltwon_arr(2:kk);

    oo_infu = log(ratio_infu)/log(2);
    oo_twou = log(ratio_twou)/log(2);
    oo_infn = log(ratio_infn)/log(2);
    oo_twon = log(ratio_twon)/log(2);
    
    linf_arr = [linfn_arr; linfu_arr];
    ltwo_arr = [ltwon_arr; ltwou_arr];
    ratio_inf_arr = [ratio_infn; ratio_infu];
    ratio_two_arr = [ratio_twon; ratio_twou];
    oo_inf_arr = [oo_infn; oo_infu];
    oo_two_arr = [oo_twon; oo_twou];

    if save_file

        mms_errors.ratio_twon = ratio_twon;
        mms_errors.ratio_infn = ratio_infn;
        mms_errors.ratio_twou = ratio_twou;
        mms_errors.ratio_infu = ratio_infu;
        mms_errors.oo_twon = oo_twon;
        mms_errors.oo_infn = oo_infn;
        mms_errors.oo_twou = oo_twou;
        mms_errors.oo_infu = oo_infu;
        mms_errors.sparsefill = sparsefill;

        filename = strcat('/Volumes/DATA/matlab/transport_verification/errors_',str,'.mat');
        save(filename,'-struct','mms_errors');

    end
    %%

    if mms_plots

        x0 = 0;
        y0 = 0;
        width = 1000;
        height = 600;

        figure(3)
        set(gcf,'Position',[x0 y0 width height],'Color','w')

        subplot(3,1,1)
        plot(nxax,ex_soln,'k','Linewidth',1.5)
        set(gca,'Fontsize',25)
        xticks([])
        xlim([min(nxax) max(nxax)])
        ylabel('$n$ (m$^{-3}$)','Interpreter','latex')
        legend('$n = n_0 + n_z\sin(k_nz^2)$','Interpreter','latex','location','northwest')

        subplot(3,1,2)
        plot(vxax,ex_solu,'k','Linewidth',1.5)
        xlim([min(nxax) max(nxax)])
        xticks([])
        set(gca,'Fontsize',25)
        ylabel('$v$ (ms$^{-1}$)','Interpreter','latex')
        legend('$v = v_0 + v_z\sin(k_vz^2)$','Interpreter','latex','location','southwest')

        subplot(3,1,3)
        plot(nxax(1:3:npts),zeros(1,length(nxax(1:3:npts))),'xk','Linewidth',1.5)
        xlim([min(nxax) max(nxax)])
        xlabel('Position (m)')
        set(gca,'Fontsize',25)
        yticks([])

        if SS

            x0 = 0;
            y0 = 0;
            width = 1000;
            height = 450;

            figure(1)
            set(gcf,'Position',[x0 y0 width height],'Color','w')
            % yyaxis left
            loglog(npts_arr,0.1./npts_arr,'.-k','Markersize',15,'Linewidth',1.2)
            hold on
            % yyaxis right
            loglog(npts_arr,ltwou_arr,'-*b','Markersize',12,'Linewidth',1.2)
            % hold on
            loglog(npts_arr,linfu_arr,'-xr','Markersize',12,'Linewidth',1.2)
            % set(gca,'yscale','log')
            xlim([min(npts_arr) max(npts_arr)])
            ylabel('Error')
            xlabel('NPTS')
            % text(npts_arr,1.0./npts_arr,{'2048','1024','512','256','128'})
            text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-40,1.0/npts_arr(3)-7.17e-3,...
                '$\mathcal{O}(\Delta z)=1$','Interpreter','latex',...
                'BackgroundColor','w','Rotation',-10)
            set(gca,'Linewidth',1.0)
            text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-49,ltwou_arr(2)+2.3e-3,...
                '$L_2$ slope $\approx$ 1','Interpreter','latex',...
                'BackgroundColor','w','Rotation',-10,'Color','b')
            text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-50,linfu_arr(2)+5.0e-3,...
                '$L_{\infty}$ slope $\approx$ 1','Interpreter','latex',...
                'BackgroundColor','w','Rotation',-10,'Color','r')
            % xticks([])
            xticks([fliplr(npts_arr)])
            set(gca,'XMinorTick','off')

        %     export_fig('/Volumes/DATA/thesis/RFT/figs/MMS_SS_mom_dz_coupled.png',...
        %         '-r300')



            figure(2)
            set(gcf,'Position',[x0 y0 width height],'Color','w')
            % yyaxis left
            loglog(npts_arr,1.0./npts_arr,'.-k','Markersize',15,'Linewidth',1.2)
            hold on
            % yyaxis right
            loglog(npts_arr,ltwon_arr,'-*b','Markersize',12,'Linewidth',1.2)
            % hold on
            loglog(npts_arr,linfn_arr,'-xr','Markersize',12,'Linewidth',1.2)
            % set(gca,'yscale','log')
            xlim([min(npts_arr) max(npts_arr)])
            ylabel('Error')
            xlabel('NPTS')
            % text(npts_arr,1.0./npts_arr,{'2048','1024','512','256','128'})
            text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-35,1.0/npts_arr(3)-1.5e-3,...
                '$\mathcal{O}(\Delta z)=1$','Interpreter','latex',...
                'BackgroundColor','w','Rotation',-18)
            set(gca,'Linewidth',1.0)
            text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-50,ltwon_arr(2)+0.6e-3,...
                '$L_2$ slope $\approx$ 1','Interpreter','latex',...
                'BackgroundColor','w','Rotation',-16,'Color','b')
            text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-48,linfn_arr(2)+1.e-3,...
                '$L_{\infty}$ slope $\approx$ 1','Interpreter','latex',...
                'BackgroundColor','w','Rotation',-16,'Color','r')
            % xticks([])
            xticks([fliplr(npts_arr)])
            set(gca,'XMinorTick','off')

        elseif TD 

            x0 = 0;
            y0 = 0;
            width = 1000;
            height = 450;

            figure(1)
            set(gcf,'Position',[x0 y0 width height],'Color','w')
            % yyaxis left
            loglog(dt_arr,dt_arr,'.-k','Markersize',15,'Linewidth',1.2)
            hold on
            % yyaxis right
            loglog(dt_arr,ltwou_arr,'-*b','Markersize',12,'Linewidth',1.2)
            % hold on
            loglog(dt_arr,linfu_arr,'-xr','Markersize',12,'Linewidth',1.2)
            % set(gca,'yscale','log')
            xlim([min(dt_arr) max(dt_arr)])
            ylabel('Error')
            xlabel('$\Delta t$','Interpreter','latex')
            % text(npts_arr,1.0./npts_arr,{'2048','1024','512','256','128'})
            text((dt_arr(3)-(dt_arr(3)-dt_arr(2))/2)-1.8e-6,dt_arr(3)-8.0e-7,...
                '$\mathcal{O}(\Delta t)=1$','Interpreter','latex',...
                'BackgroundColor','w','Rotation',6)
            set(gca,'Linewidth',1.0)
            text((dt_arr(3)-(dt_arr(3)-dt_arr(2))/2)-1.87e-6,ltwou_arr(2)-2.1e-1,...
                '$L_2$ slope $\approx$ 1','Interpreter','latex',...
                'BackgroundColor','w','Rotation',4,'Color','b')
            text((dt_arr(3)-(dt_arr(3)-dt_arr(2))/2)-1.9e-6,linfu_arr(2)-7.0e-2,...
                '$L_{\infty}$ slope $\approx$ 1','Interpreter','latex',...
                'BackgroundColor','w','Rotation',4,'Color','r')
            xticks([fliplr(dt_arr)])
            xticklabels([{'5\times10^{-7}','1\times10^{-6}','2\times10^{-6}',...
                '4\times10^{-6}','8\times10^{-6}'},'Interpreter','latex'])
        %     set(gca,'XMinorTick','off')

        %     export_fig('/Volumes/DATA/thesis/RFT/figs/MMS_TD_mom_coupled.png',...
        %         '-r300')



            figure(2)
            set(gcf,'Position',[x0 y0 width height],'Color','w')
            % yyaxis left
            loglog(dt_arr,dt_arr,'.-k','Markersize',15,'Linewidth',1.2)
            hold on
            % yyaxis right
            loglog(dt_arr,ltwon_arr,'-*b','Markersize',12,'Linewidth',1.2)
            % hold on
            loglog(dt_arr,linfn_arr,'-xr','Markersize',12,'Linewidth',1.2)
            % set(gca,'yscale','log')
            xlim([min(dt_arr) max(dt_arr)])
            ylabel('Error')
            xlabel('$\Delta t$ (s)','Interpreter','latex')
            % text(npts_arr,1.0./npts_arr,{'2048','1024','512','256','128'})
            text((dt_arr(3)-(dt_arr(3)-dt_arr(2))/2)-5.0e-7,dt_arr(3)+5.0e-7,...
                '$\mathcal{O}(\Delta t)=1$','Interpreter','latex',...
                'BackgroundColor','w','Rotation',5)
            set(gca,'Linewidth',1.0)
            text((dt_arr(3)-(dt_arr(3)-dt_arr(2))/2)-7.2e-7,ltwon_arr(2)-1.5e-1,...
                '$L_2$ slope $\approx$ 1','Interpreter','latex',...
                'BackgroundColor','w','Rotation',5,'Color','b')
            text((dt_arr(3)-(dt_arr(3)-dt_arr(2))/2)-8.0e-7,linfn_arr(2)+8e-2,...
                '$L_{\infty}$ slope $\approx$ 1','Interpreter','latex',...
                'BackgroundColor','w','Rotation',4,'Color','r')
            % xticks([])
            xticks([fliplr(dt_arr)])
            xticklabels([{'5\times10^{-7}','1\times10^{-6}','2\times10^{-6}',...
                '4\times10^{-6}','8\times10^{-6}'},'Interpreter','latex'])
        %     set(gca,'XMinorTick','off')    
        end

    end

end



