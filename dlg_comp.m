%-----------------------------------------%
% parameters for comparison with          %
% dlg code                                %
% rlbarnett c3149416, 151217              %
%-----------------------------------------%

%--
% driver terms
% freq = 42.0e6;
freq = 5.9973e+09;
om = 2*pi*freq;
source_width = 0.01;
source_loc = 1.45;
dampFac = 10.0;

%--
% define wavenumbers ky and kz (/m); use values given in van eester section IV?? No
% others mentioned
ky = 10.0 + 10.0i;
kz = 0.0;
k0 = om/const.c0;

%--
% "common local derivatives for N0, v||^2, static potential and
% ponderomotive potential" 
lamby = 0.0;
lambz = 0.0;

%--
% spatial domain
npts = 256;
xmin = 1.2;
xmax = 1.7;
dx = (xmax - xmin)/(npts - 1);
% npts = ((xmax - xmin)/dx);
xax = linspace(xmin, xmax, npts);
% xax = xmin:dx:xmax;
% dx = (xmax - xmin)/(npts - 1);

%--
% magnetic field (tesla)
R0 = 1.32;
% B0 = 2.6*R0./xax;
B0 = 0.0;

%--
% temperature
T_ev = 15.0;

%--
% background density -- set to zero for vacuum case
Nmax = 50.0e19;
% Nmin = 50.0e16;
% N0 = Nmax*ones(1,npts);
N0 = 0.0;
% m = (Nmax - Nmin) ./ (xmax - xmin);
% N0 = 10.^(m*xax + Nmin);

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

np_bound = floor(0.2*npts);
ax = linspace(0,pi,np_bound);
damp0 = (cos(ax)+1)/2;
damp = ones(1,npts);
damp(1:np_bound) = damp(1:np_bound) + dampFac*i*damp0;
% damp(1:end) = damp(end-(np_bound-1):end) + dampFac*i*fliplr(damp0);

me = me .* damp;

%--
% densities
N0e = 1.0*N0;
N0d = 0.95*N0;
N0h = 0.05*N0;
N0i = N0h + N0d;






