%-----------------------------------------%
% parameters for comparison with          %
% dve 2015 paper                          %
% rlbarnett c3149416, 201217              %
%-----------------------------------------%

%--
% driver terms
freq = 51.0e6;
source_width = 0.001;
source_loc = 0.19;
dampFac = 10.0;

%--
% define wavenumbers ky and kz (/m); use values given in van eester section IV?? No
% others mentioned
ky = 0.0;
kz = 0.0;

%--
% "common local derivatives for N0, v||^2, static potential and
% ponderomotive potential" 
lamby = 0.0;
lambz = 0.0;

%--
% spatial domain
npts = 256;
xmin = 0.0;
xmax = 0.2;
dx = (xmax - xmin)/(npts - 1);
% npts = ((xmax - xmin)/dx);
xax = linspace(xmin, xmax, npts);
% xax = xmin:dx:xmax;
% dx = (xmax - xmin)/(npts - 1);

%--
% magnetic field (tesla)
B0 = 2.6*ones(1,npts);

%--
% temperature
T_ev = 15.0;

%--
% background density -- set to zero for vacuum case
Nmax = 17;
Nmin = 12;
m = (Nmax - Nmin) ./ (xmax - xmin);
N0 = 10.^(m*xax + Nmin);
% N0 = Nmax*ones(1,npts);

%--
% initialise perturbed density as zero
N1 = zeros(1,npts);
N1e = N1;
N1i = N1;

%--
% initialise perturbed velocity as zero
v1 = zeros(1,npts);
v1e = v1;
v1i = v1;

%--
% electron mass
me = 9.11e-31*ones(1,npts);

% DLG - since I don't have the license for "makedist"
% I fixed your cos ramping function :)

np_bound = floor(0.07*npts);
ax = linspace(0,pi,np_bound);
damp0 = (cos(ax)+1)/2;
damp = ones(1,npts);
% damp(1:np_bound) = damp(1:np_bound) + dampFac*1i*damp0;
damp(end-(np_bound-1):end) = damp(end-(np_bound-1):end) + dampFac*1i*fliplr(damp0);

me = me .* damp;

%--
% densities
N0e = 1.0*N0;
N0d = 0.95*N0;
N0h = 0.05*N0;
N0i = N0h + N0d;

