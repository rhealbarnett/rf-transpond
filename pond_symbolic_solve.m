syms ax ay az ux uy uz OMz om q m Ex Ey Ez

q = sym('q', 'real');
om = sym('om', 'real');
m = sym('m', 'real');

eq1 = -1i*om*ux - uy*OMz - ax == 0;
eq2 = om*uy + 1i*ux*OMz - 1i*ay == 0;

S = solve(eq1,eq2);

ux = S.ux;
uy = S.uy;
uz = 1i*az/om;

ax = q*Ex/m;
ay = q*Ey/m;
az = q*Ez/m;

ux = subs(ux); uy = subs(uy); uz = subs(uz);

syms Dx DDx D S P kz k0 h k

Roty = [floor(cos(pi/2)), 0 sin(pi/2);0,1,0;sin(pi/2),0,floor(cos(pi/2))];
K = [kz^2, 0, 1i*kz*Dx;
    0, kz^2 - DDx, 0;
    1i*kz*Dx, 0, -DDx];

Kdiff = [0,0,-1i*k/(2*h);
         0,-1/h^2,0;
         -1i*k/(2*h),0,-1.0/h^2];
    

