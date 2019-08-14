syms ax ay az ux uy uz OMz om q m Ex Ey Ez

q = sym('q', 'real');
om = sym('om', 'real');
m = sym('m', 'real');

eq1 = -1i*om*ux - uy*OMz - ax == 0;
eq2 = -1i*om*uy + ux*OMz - ay == 0;

S = solve(eq1,eq2);

ux = S.ux;
uy = S.uy;
uz = 1i*az/om;

ax = q*Ex/m;
ay = q*Ey/m;
az = q*Ez/m;

ux = subs(ux); uy = subs(uy); uz = subs(uz);

syms ky S D k0

eq1 = (ky^2 - k0^2)*Ex + (1i*ky + 1i*D)*Ey == 0;

S = solve(eq1,Ey);