%------------------------------------------%
% calculate pressure term                  %
% -vt^2*del*N0/N0                          %
% rlbarnett c3149416 141117                %
%------------------------------------------%

gradN0ey = lamby*N0e;
gradN0ez = lambz*N0e;

gradN0iy = lamby*N0i;
gradN0iz = lambz*N0i;

% gradN0e = zeros(npts,1);
% gradN0i = zeros(npts,1);
% 
% for ii=2:(npts-1)
%     
%     gradN0e(ii) = (N0e(ii-1) - N0e(ii+1))/(2.0*dx);
%     gradN0i(ii) = (N0i(ii-1) - N0i(ii+1))/(2.0*dx);
%     
% end
% 
% gradN0e(1) = (-3.0*N0e(1) + 4.0*N0e(2) - N0e(3))/(2.0*dx);
% gradN0e(npts) = (3.0*N0e(npts) - 4.0*N0e(npts-1) + N0e(npts-2))/(2.0*dx);
% 
% gradN0i(1) = (-3.0*N0i(1) + 4.0*N0i(2) - N0i(3))/(2.0*dx);
% gradN0i(npts) = (3.0*N0i(npts) - 4.0*N0i(npts-1) + N0i(npts-2))/(2.0*dx);

gradN0ex = gradient(N0e, dx);
gradN0ix = gradient(N0i, dx);

pressex = -vt^2*gradN0ex/N0e;
pressey = -vt^2*gradN0ey/N0e;
pressez = -vt^2*gradN0ez/N0e;

pressix = -vt^2*gradN0ix/N0i;
pressiy = -vt^2*gradN0iy/N0i;
pressiz = -vt^2*gradN0iz/N0i;

presse = [pressex, pressey, pressez];
pressi = [pressix, pressiy, pressiz];
