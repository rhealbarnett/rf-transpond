%------------------------------%
% plots for thesis 
% mostly white space format!
%------------------------------%

function ans = thesis_fig(ax,variable,labely,labelx,line_width,str)

    x0 = 0;
    y0 = 0;
    width = 900;
    height = 400;

    figure(1)
    set(gcf,'Position',[x0 y0 width height],'Color','w')
    plot(ax,variable,'k','Linewidth',line_width)
    xlim([min(ax) max(ax)])
    xlabel(labelx,'Interpreter','latex')
    ylabel(labely,'Interpreter','latex')
    
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