%--------------------------------------------------------------------------------------------------------------%
% Input file for transport 1.5D                                                                            %
% nabla = (kapx, kapy, dz)
% rlbarnett c3149416 190723     
%--------------------------------------------------------------------------------------------------------------%
%%

%-------
% call file containing physical constants. %
%-------
const = constants();
mp = const.mp;
e = const.e;

%-------
% define parameters. %
%-------
Te = 5.0;
Ti = 0.5;
cs = sqrt((Te + Ti)*e/m);

eta_para = 1.0;

%-- 
% D_perp will be used to calculate eta_perp once the density (n) 
% has been defined. 
D_perp = 1.0;

%-- 
% Setting the perp derivatives as zero to begin, for comparison with the 
% parallel only case. 
kapx = 0.0;
kapy = 0.0;


