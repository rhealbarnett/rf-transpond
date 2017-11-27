%------------------------------------------%
% dydx function                            %
% coefficients A(x), B(x), C(x)            %
% rlbarnett c3149416 211117                %
%------------------------------------------%

%--
% formulate differential equation
dydx = @(x,y) -(B(x).*y + C(x))/A(x);

%--
% find solution
sol = bvp5c(dydx,bound,solinit);








