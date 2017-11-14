%------------------------------------------%
% calculate pressure term                  %
% -vt^2*del*N0/N0                          %
% rlbarnett c3149416 141117                %
%------------------------------------------%

gradNe = zeros(npts,1);
gradNi = zeros(npts,1);

for ii=2:(npts-1)
    
    gradNe(ii) = (Ne(ii-1) - Ne(ii+1))/(2.0*dx);
    gradNi(ii) = (Ni(ii-1) - Ni(ii+1))/(2.0*dx);
    
end

gradNe(1) = (-3.0*Ne(1) + 4.0*Ne(2) - Ne(3))/(2.0*dx);
gradNe(npts) = (3.0*Ne(npts) - 4.0*Ne(npts-1) + Ne(npts-2))/(2.0*dx);

gradNi(1) = (-3.0*Ni(1) + 4.0*Ni(2) - Ni(3))/(2.0*dx);
gradNi(npts) = (3.0*Ni(npts) - 4.0*Ni(npts-1) + Ni(npts-2))/(2.0*dx);