%------------------------------------------%
% ODE test solve                           %
% bvp4c for boundary value problems        %
% rlbarnett c3149416 171117                %
%------------------------------------------%

xax = linspace(1,2,10);
bound = 1.0;

A_test = 2.0.*xax.^2;
B_test = 3.0.*xax;
C_test = 1.0./xax;

A = @(x) interp1(xax,A_test,x);
B = @(x) interp1(xax,B_test,x);
C = @(x) interp1(xax,C_test,x);

dydx = @(x,y) -(B(x)*y + C(x))/A(x);
bounds = @(ya,yb) ya - bound;

solinit = bvpinit(xax,[1]);

sol = bvp4c(dydx,bounds,solinit);


% function dydx = ode(x,y)
%     
%     dydx = y(1);
% %     dydx = [y(2); -abs(y(1))];
%     
% end
% 
% function res = bcs(ya,yb,bound)
% 
%     res = [ya(1) - bound];%; yb(1) - 0.3679];
% %     res = [ya(1); yb(1) + 2];
%     
% end


