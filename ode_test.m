%------------------------------------------%
% ODE test solve                           %
% bvp4c for boundary value problems        %
% rlbarnett c3149416 171117                %
%------------------------------------------%

solinit = bvpinit([0 1],[0]);

bound = 1.0;

sol = bvp4c(@ode,@(ya,yb)bcs(ya,yb,bound),solinit);


function dydx = ode(x,y)
    
    dydx = y(1);
%     dydx = [y(2); -abs(y(1))];
    
end

function res = bcs(ya,yb,bound)

    res = [ya(1) - bound];%; yb(1) - 0.3679];
%     res = [ya(1); yb(1) + 2];
    
end


