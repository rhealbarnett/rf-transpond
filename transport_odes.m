%----------------------------------------%
% solve coupled set of transport ODEs    %
% dvz/dt + (1/mn)(Te+Ti)dn/dz = 0        %
% dn/dt + (n)dvz/dz = 0                  %
% rlbarnett c3149416 050318              %
%----------------------------------------%

syms v(t,z) n(t,z)

m = 1.67e-27;
Te = 5.0;
Ti = Te;

eq1 = diff(v,t) + (1.0/(m*n(t,z)))*(Te + Ti)*diff(n,z) == 0.0;
eq2 = diff(n,t) + n(t,z)*diff(v,z) == 0.0;
eqns = [eq1, eq2];

[sol1, sol2] = dsolve(eqns);