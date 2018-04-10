%-----------------------------------------%
% mms for comparison with                 %
% dlg code                                %
% rlbarnett c3149416, 090418              %
%-----------------------------------------%

%--
% constants
mu0 = 4.0*pi*1.0e-7;
eps0 = 8.85e-12;
c0 = 1.0/sqrt(eps0*mu0);

%--
% driver freq
freq = 13.0e6;
om = 2.0*pi*freq;
k0 = om/c0;
wavel0 = (2*pi)/k0;

%
npts = 64;
xmin = -1.0;
xmax = 1.0;

%--
% electron constants
e = -1.6022e-19;
me = 9.11e-31*ones(1,npts);

%--
% ion constants (95% D, 5% H)
qd = abs(e);
mp = 1.67e-27;
md = 2.0*mp*ones(1,npts);

qh = abs(e);
mh = mp*ones(1,npts);

B0 = 0.0;
N0 = 0.0;

ky = 0.0;
kz = 0.0;

%-- 
% cyclotron frequencies
om_ce = e*B0./me;
om_cd = qd*B0./md;
om_ch = qh*B0./mh;

%--
% plasma frequencies
om_pe = sqrt(N0*e^2./(me*eps0));
om_pd = sqrt(N0*qd^2./(md*eps0));
om_ph = sqrt(N0*qh^2./(mh*eps0));

%--
% rotation matrix
alpha = 0.0;
beta = 0.0;

r11 = cos(beta)*cos(alpha);
r12 = cos(beta)*sin(alpha);
r13 = -sin(beta);
r21 = -sin(alpha);
r22 = cos(alpha);
r23 = 0.0;
r31 = sin(beta)*cos(alpha);
r32 = sin(beta)*sin(alpha);
r33 = cos(beta);

rot = [[r11, r12, r13]
     [r21, r22, r23]
     [r31, r32, r33]];
 
e_para = rot(3,:);
Bvec = B0.*transpose(repmat(e_para,npts,1));

dielec_tens;