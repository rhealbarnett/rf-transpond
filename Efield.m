%-----------------------------------------------------%
% analytic E field for PF in transport code
% based loosely on comsol LAPD results for E||
% rlbarnett c3149416 190304
%-----------------------------------------------------%



sigma = sqrt(0.01);
mu = (xmax-abs(xmin))/2;
efield = @(x) (1.0/sqrt(2.0*pi*sigma^2))*exp(-(x - mu).^2/(2.0*sigma^2));

fwhm = 2.0*sqrt(2.0*log(2))*sigma;