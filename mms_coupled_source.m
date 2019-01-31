syms n0 u0 nx ux knx kux x lam om t 

n = n0 + nx*sin(knx*x^2 + om*t);

u = u0 + ux*cos(kux*x^2 + om*t);

nu = n*u;

dndt = diff(n,t);
dnudx = diff(nu,x);
dudt = diff(u,t);
dudx = diff(u,x);
d2udx = diff(u,x,2);

