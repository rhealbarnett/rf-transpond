%------------------------------%
% plots for thesis 
% mostly white space format!
%------------------------------%

function ans = thesis_fig(num,ax,variable,labely,labelx,line_width,colour,str)

    x0 = 0;
    y0 = 0;
    width = 1000;
    height = 500;
    
    col_purp = [0.28 0.06 0.47];
    col_pink = [0.78 0.24 0.3];
    col_oran = [0.94 0.44 0.13];

    figure(num)
    set(gcf,'Position',[x0 y0 width height],'Color','w')
    plot(ax,variable,'color',colour,'Linewidth',line_width)
    xlim([min(ax) max(ax)])
    xlabel(labelx,'Interpreter','latex')
    ylabel(labely,'Interpreter','latex')
%     hold on
    
%     dim = [.47 .7 .25 .15];
%     annotation('textbox',dim,'String',str,'Interpreter','latex',...
%         'Edgecolor','none')
%     dim = [.73 .7 .25 .15];
%     annotation('textbox',dim,'String','$v = c_s$ ms$^{-1}$','Interpreter','latex',...
%         'Edgecolor','none')
%     dim = [.18 .15 .25 .15];
%     annotation('textbox',dim,'String','$v = -c_s$ ms$^{-1}$','Interpreter','latex',...
%         'Edgecolor','none')

end