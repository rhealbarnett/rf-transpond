syms ax ay az ux uy uz OMz om q m Ex Ey Ez Z e n PI

eq1 = -1i*om*ux - uy*OMz - ax == 0;
eq2 = om*uy + 1i*ux*OMz - 1i*ay == 0;

S = solve(eq1,eq2);

ux = S.ux;
uy = S.uy;

ax = q*Ex/m;
ay = q*Ey/m;
az = q*Ex/m;

% ux = simplify(subs(ux));
% uy = simplify(subs(uy));

coeff = (4*PI*1i/om)*Z*e*n;

term1 = diff(conj(Ex(x)),x)*ux + diff(conj(Ey(x)),x)*uy + diff(conj(Ez(x)),x)*uz;
term2 = diff(ux(x)*conj(Ex),x);


% syms nj njp1 dx vj vjm1
% 
% check = (1.0/nj)*((njp1 - nj)/dx)*((vj-vjm1)/dx);
% 
% syms dx dy dz Ex Ey Ez Dx Dy Dz ux uy uz
% 
% nabla = [dx; dy; dz];
% E = [Ex; Ey; Ez];
% D = [Dx; Dy; Dz];
% u = [ux; uy; uz];
% 
% term1 = (nabla*conj(transpose(E)))*D;
% term2 = transpose(nabla)*(D*conj(transpose(E)));
% term3 = transpose(nabla)*(u*transpose(u));