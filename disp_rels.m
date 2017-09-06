%-----------------------------------------%
% Dispersion relation code                %
% Adapted from python for comparison      %
% Determine wtf the eigenvecs mean??      %
% rlbarnett c3149416, 060917              %
%-----------------------------------------%

%-- magnetic field
B0 = 1.0;

%-- 'normalised' mass? And mi/me = 5.0
me = 1.0;
mi = 5.0*me;

%-- 'normalised' constants? 
eps0 = 1.0;
mu0 = 1.0;
c0 = 1.0/sqrt(mu0*eps0);
e = -1.0;
qi = abs(e);

%-- density
alpha = 3.0;
N = alpha*eps0*(B0^2)/me;

%-- plasma frequencies
om_pe = sqrt(N*e^2/(eps0*me));
om_pi = sqrt(N*qi^2/(eps0*mi));

%-- cyclotron frequencies
om_ce = e*B0/me;
om_ci = qi*B0/mi;

%-- dielectric tensor elements
s = 1.0 - om_pe^2/(om^2 - om_ce^2) - om_pi^2/(om^2 - om_ci^2);
d = om_ce*om_pe^2/(om*(om^2 - om_ce^2)) + om_ci*om_pi^2/(om*(om^2 - om_ci^2));
p = 1.0 - om_pe^2/om^2 - om_pi^2/om^2;

s = s + 1.0*j;
d = d + 1.0*j;
p = p + 1.0*j;

r = s + d;
l = s - d;

cpdt = [[s, -j*d, 0.0],
        [j*d, s, 0.0],
        [0.0, 0.0, p]];
    
theta = pi/2.0;

a = s*sin(theta)^2 + p*cos(theta)^2;
b = r*l*sin(theta)^2 + p*s*(1.0 + cos(theta)^2);
c = p*r*l

f = sqrt(b^2 - 4.0*a*c)

ns_pos = (b + f)/(2.0*a)
ns_neg = (b + f)/(2.0*a)




