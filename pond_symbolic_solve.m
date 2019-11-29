syms ax ay az ux uy uz OMz om q m Ex Ey Ez

q = sym('q', 'real');
om = sym('om', 'real');
m = sym('m', 'real');

eq1 = -1i*om*ux + uy*OMz - ax == 0;
eq2 = om*uy - 1i*ux*OMz - 1i*ay == 0;

S = solve(eq1,eq2);

ux = S.ux;
uy = S.uy;
uz = 1i*az/om;

ax = q*Ex/m;
ay = q*Ey/m;
az = 0;

ux = subs(ux); uy = subs(uy); uz = subs(uz);

clear all
syms ux(z) uy(z) uz(z) Jx(z) Jy(z) Jz(z) Ex(z) Ey(z) Ez(z)...
    q om m n az(z) OMz ax(z) ay(z) x y z

E = [Ex(z), Ey(z), 0];
J = [Jx(z), Jy(z), 0];
u = [ux(z), uy(z), 0];

JE = mtimes(transpose(J),conj(E));
uu = mtimes(transpose(u),conj(u));
divJE = sym('divJE',[3,1]);
divuu = sym('divuu',[3,1]);
gradE = sym('gradE',[3,3]);
gradEdotJ = sym('gradEdotJ',[3,1]);
vec = [x,y,z];

for ii=1:3
    
    divJE(ii,1) = diff(JE(1,ii),x) + diff(JE(2,ii),y) + diff(JE(3,ii),z);
    divuu(ii,1) = diff(uu(1,ii),x) + diff(uu(2,ii),y) + diff(uu(3,ii),z);
    
    for jj=1:3
        
        gradE(jj,ii) = diff(conj(E(1,ii)),vec(jj));
        
    end
   
end

for ii=1:3
    
    gradEdotJ(ii,1) = J(1)*gradE(ii,1) + J(2)*gradE(ii,2) + J(3)*gradE(ii,3);
    
end

pf_vec = (1i/om)*(gradEdotJ - divJE) - n*m*divuu;

pf_vec = subs(pf_vec,J,n*q*u);
pf_vec = subs(pf_vec,[ux(z),uy(z)],[(1i*om*ax(z)-OMz*ay(z))/(om^2 - OMz^2),...
    (OMz*ax(z) + 1i*om*ay(z))/(om^2 - OMz^2)]);
% pf_vec = subs(pf_vec,uz(z),(1i/om)*az(z));
% pf_vec = subs(pf_vec,az(z),q*Ez(z)/m);
pf_vec = subs(pf_vec,[ax(z), ay(z)],[q*Ex(z)/m, q*Ey(z)/m]);

assume(q,'real')
assume(om,'real')
assume(m,'real')
assume(n,'real')
assume(OMz,'real')

PF = (1.0/2.0)*real(pf_vec);

assume(q,'real')
assume(om,'real')
assume(m,'real')
assume(n,'real')
assume(OMz,'real')


 -(imag((n*q*conj(diff(Ex(z), z))*((OMz*q*Ey(z))/m - (om*q*Ex(z)*1i)/m))/(OMz^2 - om^2)) -...
     imag((n*q*conj(diff(Ey(z), z))*((OMz*q*Ex(z))/m + (om*q*Ey(z)*1i)/m))/(OMz^2 - om^2)))/(2*om);

