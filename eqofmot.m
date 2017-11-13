%------------------------------------------%
% solve eqn of motion                      %
% adapted from python code 130417          %
% rlbarnett c3149416 301017                %
%------------------------------------------%

%--
% number of dimensions (spatial + velocity)
pos_dim = 3;
n = 2;

%--
% initial velocity 
v0 = 0.1*vt;

vx = v0;
vy = v0;
vz = v0;

v = [vx, vy, vz];

%--
% initial position
xx = 0.0;
xy = 0.0;
xz = 0.0;

x = [xx, xy, xz];

%--
% initialise arrays for position, velocity and acceleration
boX_arr = zeros(nmax, pos_dim, n);
a_arr = zeros(nmax, pos_dim);

%--
% BORIS PUSH ALGORITHM
% magnetic field rotation
OM = (e / me)*(dt / 2.0)*Bvec;
OM_mag_sq = OM(1)^2.0 + OM(2)^2.0 + OM(3)^2.0;
s = 2.0*OM / (1.0 + OM_mag_sq);

for ii=1:nmax
    
    boX_arr(ii,:,1) = x;
    boX_arr(ii,:,2) = v;
    
    ax = 1.0e6*(8.0*x(1)^3.0);
    ay = 1.0e6*(20.0*x(2)^3.0);
    az = 1.0e6*(40.0*x(3)^3.0);
    
    a = [ax, ay, az];
    a_arr(ii,:) = a;
    
    aE = (dt / 2.0)*a;
    
    vm = v + aE;
    vpr = vm + cross(vm, OM);
    vpl = vm + cross(vpr, s);
    v = vpl + aE;
    
    x = x + v*dt;
    
end
    




