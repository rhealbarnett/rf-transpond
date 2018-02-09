%-----------------------------------------%
% momentum & continuity solver: pdepe     %
% rlbarnett c3149416                      %
% 090218                                  %
%-----------------------------------------%

function pde_cont_mom

% mi = 1.67e-27;
% Tev = 10.0;
% cs = sqrt(Tev/mi);

m = 0;
x = linspace(0,20,20);
t = linspace(0,20,20);

sol = pdepe(m,@pdezpe,@pdezic,@pdezbc,x,t);

function [c,f,s] = pdezpe(x,t,u,DuDx)
mi = 1.67e-27;
Tev = 10.0;
c = [1; 1];
f = [(1.0/2.0)*u(1)^2; -u(2)];
s = [-1.0/(mi*u(2))*(DuDx(2))*Tev; (-u(2)/u(1))*DuDx(1)];

function u0 = pdezic(x)
u0 = [0; 1.0e19];

function [pl,ql,pr,qr] = pdezbc(xl,ul,xr,ur,t)
mi = 1.67e-27;
Tev = 10.0;
e = 1.6022e-19;
cs = sqrt(Tev*e/mi);
pl = [ul(1); ul(2) - 1.0e19];
ql = [0; 0];
pr = [ul(1) - cs; 0];
qr = [0; 1];



