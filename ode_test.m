%------------------------------------------%
% ODE test solve                           %
% bvp4c for boundary value problems        %
% rlbarnett c3149416 171117                %
%------------------------------------------%

solinit = bvpinit([0 1],[0]);

sol = bvp4c(@ode,@bcs,solinit);


function dydx = ode(x,y)
    
    dydx = -y(1);
    dydx
%     dydx = [y(2); -abs(y(1))];
    
end

function res = bcs(ya,yb)

    res = [ya(1) - 1.0];%; yb(1) - 0.3679];
%     res = [ya(1); yb(1) + 2];
    
end


