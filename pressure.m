%------------------------------------------%
% calculate pressure term                  %
% -vt^2*del*N0/N0                          %
% rlbarnett c3149416 141117                %
%------------------------------------------%

gradN0e = zeros(npts,1);
gradN0i = zeros(npts,1);

for ii=2:(npts-1)
    
    gradN0e(ii) = (N0e(ii-1) - N0e(ii+1))/(2.0*dx);
    gradN0i(ii) = (N0i(ii-1) - N0i(ii+1))/(2.0*dx);
    
end

gradN0e(1) = (-3.0*N0e(1) + 4.0*N0e(2) - N0e(3))/(2.0*dx);
gradN0e(npts) = (3.0*N0e(npts) - 4.0*N0e(npts-1) + N0e(npts-2))/(2.0*dx);

gradN0i(1) = (-3.0*N0i(1) + 4.0*N0i(2) - N0i(3))/(2.0*dx);
gradN0i(npts) = (3.0*N0i(npts) - 4.0*N0i(npts-1) + N0i(npts-2))/(2.0*dx);