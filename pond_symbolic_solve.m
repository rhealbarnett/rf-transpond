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
syms ux(x,y,z) uy(x,y,z) uz(x,y,z) Jx(x,y,z) Jy(x,y,z) Jz(x,y,z) Ex(x,y,z)...
    Ey(x,y,z) Ez(x,y,z) q om m n az(x,y,z) OMz ax(x,y,z) ay(x,y,z) x y z

E = [Ex(x,y,z), Ey(x,y,z), Ez(x,y,z)];
J = [Jx(x,y,z), Jy(x,y,z), Jz(x,y,z)];
u = [ux(x,y,z), uy(x,y,z), uz(x,y,z)];

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
% pf_vec = subs(pf_vec,[ux(z),uy(z)],[(1i*om*ax(z)-OMz*ay(z))/(om^2 - OMz^2),...
%     (OMz*ax(z) + 1i*om*ay(z))/(om^2 - OMz^2)]);
% pf_vec = subs(pf_vec,uz(z),(1i/om)*az(z));
% pf_vec = subs(pf_vec,az(z),q*Ez(z)/m);
% pf_vec = subs(pf_vec,[ax(z),ay(z)],[q*Ex(z)/m, q*Ey(z)/m]);
pf_vec = subs(pf_vec,[ux(x,y,z),uy(x,y,z),uz(x,y,z)],[(1i*om*ax(x,y,z)-OMz*ay(x,y,z))/(om^2 - OMz^2),...
    (OMz*ax(x,y,z) + 1i*om*ay(x,y,z))/(om^2 - OMz^2), (1i/om)*az(x,y,z)]);
pf_vec = subs(pf_vec,[ax(x,y,z),ay(x,y,z),az(x,y,z)],[q*Ex(x,y,z)/m,...
    q*Ey(x,y,z)/m, q*Ez(x,y,z)/m]);

assume(q,'real')
assume(om,'real')
assume(m,'real')
assume(n,'real')
assume(OMz,'real')

PF = (1.0/2.0)*real(pf_vec);
% 
% assume(q,'real')
% assume(om,'real')
% assume(m,'real')
% assume(n,'real')
% assume(OMz,'real')


 - (((n*q*conj(Ey(x, y, z))*((OMz*q*diff(Ey(x, y, z), x))/m - (om*q*diff(Ex(x, y, z), x)*1i)/m))/(OMz^2 - om^2) - 
 (n*q*conj(Ey(x, y, z))*((OMz*q*diff(Ex(x, y, z), y))/m + (om*q*diff(Ey(x, y, z), y)*1i)/m))/(OMz^2 - om^2) + 
 (n*q*conj(diff(Ey(x, y, z), x))*((OMz*q*Ey(x, y, z))/m - (om*q*Ex(x, y, z)*1i)/m))/(OMz^2 - om^2) - 
 (n*q*conj(diff(Ex(x, y, z), y))*((OMz*q*Ey(x, y, z))/m - (om*q*Ex(x, y, z)*1i)/m))/(OMz^2 - om^2) - 
 (n*q^2*conj(diff(Ez(x, y, z), y))*Ez(x, y, z)*1i)/(m*om) + 
 (n*q^2*conj(Ey(x, y, z))*diff(Ez(x, y, z), z)*1i)/(m*om) + 
 (n*q^2*conj(diff(Ey(x, y, z), z))*Ez(x, y, z)*1i)/(m*om))*1i)/om + 
     m*n*((((q*conj(Ex(x, y, z))*conj(OMz))/m - 
     (q*conj(Ey(x, y, z))*conj(om)*1i)/m)*((OMz*q*diff(Ey(x, y, z), x))/m - 
     (om*q*diff(Ex(x, y, z), x)*1i)/m))/((OMz^2 - om^2)*(OMz^2 - om^2)) + 
     (((OMz*q*Ey(x, y, z))/m - (om*q*Ex(x, y, z)*1i)/m)*((q*conj(OMz)*conj(diff(Ex(x, y, z), x)))/m - 
     (q*conj(om)*conj(diff(Ey(x, y, z), x))*1i)/m))/((OMz^2 - om^2)*(OMz^2 - om^2)) - 
     (((q*conj(OMz)*conj(diff(Ex(x, y, z), y)))/m - 
     (q*conj(om)*conj(diff(Ey(x, y, z), y))*1i)/m)*((OMz*q*Ex(x, y, z))/m + 
     (om*q*Ey(x, y, z)*1i)/m))/((OMz^2 - om^2)*(OMz^2 - om^2)) - 
     (((q*conj(Ex(x, y, z))*conj(OMz))/m - 
     (q*conj(Ey(x, y, z))*conj(om)*1i)/m)*((OMz*q*diff(Ex(x, y, z), y))/m + 
     (om*q*diff(Ey(x, y, z), y)*1i)/m))/((OMz^2 - om^2)*(OMz^2 - om^2)) + 
     (q*((q*conj(OMz)*conj(diff(Ex(x, y, z), z)))/m - 
     (q*conj(om)*conj(diff(Ey(x, y, z), z))*1i)/m)*Ez(x, y, z)*1i)/(m*om*(OMz^2 - om^2)) + 
     (q*((q*conj(Ex(x, y, z))*conj(OMz))/m - 
     (q*conj(Ey(x, y, z))*conj(om)*1i)/m)*diff(Ez(x, y, z), z)*1i)/(m*om*(OMz^2 - om^2)))

