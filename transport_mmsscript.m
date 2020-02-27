% npts = 4096*4;
% npts = 128;
% pow = 7;
npts = 512;

transport_mms;

% dt = 8.0e-6;
dt = 1.0;
tol = 1.0e-14;
% tmin = 0;
% tmax = 8.0e-4;

for kk=1:5
    
%     npts = 2^pow;
%     nmax = round(tmax/dt);
    nmax = 1.0e5;
    save_iter = round(nmax/5);
    
    fprintf('iteration %d\n', kk)
    fprintf('Number of points %d\n', npts)
%     fprintf('Time step %d\n', dt)
%     fprintf('nmax, tmax %d, %d\n', [nmax tmax])
    
    transport_1d
    
    fprintf('Minimum dx %d\n', min(ndx))
    fprintf('Maximum dx %d\n', max(ndx))
    
    ltwon_arr(1,kk) = l_twon;
    linfn_arr(1,kk) = l_infn;
    ltwou_arr(1,kk) = l_twou;
    linfu_arr(1,kk) = l_infu;
    
%     pow = pow + 1;
    npts_arr(1,kk) = npts;

%     dt_arr(1,kk) = dt;
%     dt = dt/2.0;
    
    npts = npts / 2;
end


ratio_infu = linfu_arr(2:kk)./linfu_arr(1:kk-1);
ratio_twou = ltwou_arr(2:kk)./ltwou_arr(1:kk-1);
ratio_infn = linfn_arr(2:kk)./linfn_arr(1:kk-1);
ratio_twon = ltwon_arr(2:kk)./ltwon_arr(1:kk-1);

oo_infu = log(ratio_infu)/log(2);
oo_twou = log(ratio_twou)/log(2);
oo_infn = log(ratio_infn)/log(2);
oo_twon = log(ratio_twon)/log(2);

%%

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
    '$\mathcal{O}(\Delta x)=1$','Interpreter','latex',...
    'BackgroundColor','w','Rotation',-10)
set(gca,'Linewidth',1.0)
text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-45,ltwou_arr(2)+2.5e-3,...
    '$L_2$ slope $\approx$ 1','Interpreter','latex',...
    'BackgroundColor','w','Rotation',-10,'Color','b')
text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-45,linfu_arr(2)+5.0e-3,...
    '$L_{\infty}$ slope $\approx$ 1','Interpreter','latex',...
    'BackgroundColor','w','Rotation',-10,'Color','r')
% xticks([])
xticks([fliplr(npts_arr)])
set(gca,'XMinorTick','off')

% export_fig('/Volumes/DATA/thesis/RFT/figs/MMS_SS_mom_dx_coupled.png',...
%     '-r300')

%%

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
    '$\mathcal{O}(\Delta x)=1$','Interpreter','latex',...
    'BackgroundColor','w','Rotation',-18)
set(gca,'Linewidth',1.0)
text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-45,ltwon_arr(2)+1.0e-3,...
    '$L_2$ slope $\approx$ 1','Interpreter','latex',...
    'BackgroundColor','w','Rotation',-18,'Color','b')
text((npts_arr(3)-(npts_arr(3)-npts_arr(2))/2)-45,linfn_arr(2)+2.0e-3,...
    '$L_{\infty}$ slope $\approx$ 1','Interpreter','latex',...
    'BackgroundColor','w','Rotation',-18,'Color','r')
% xticks([])
xticks([fliplr(npts_arr)])
set(gca,'XMinorTick','off')

% export_fig('/Volumes/DATA/thesis/RFT/figs/MMS_SS_cont_dx_coupled.png',...
%     '-r300')

