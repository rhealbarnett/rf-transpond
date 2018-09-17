%%This function calculates Coloumb collisions given an ion mass mi (kg),
%%an ion charge Z,an ion temperature Ti (Please use J here), Te (Please use
%%Joules here), and electron density (please use SI m.^-3)

%%%%%%%%%%%% IF YOU DIDN'T READ ABOVE MAKE SURE YOU USE ALL S.I. UNITS (YES
%%%%%%%%%%%% THIS REQUIRES TEMPERATURE TO BE DEFINED IN Joules NOTE eV; 
%%%%%%%%%%%% BLASPHEMUS, I KNOW.)

function [vei, vee, vie, vii] = Coloumb_collisions_v2(mi,Z,Ti,Te,ne)

%%%%%Redefine constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constants;
e_const = e;
epsilon0_const = eps0;
me = me;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Define inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mi = 2*amu;
% Z = 1;
% Ti = 2*q;
% Te = 3*q;
% ne = 3e19;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mree = (me*me)./(me+me);
mrei = (me.*mi)./(me+mi);
mrii = (mi.*mi)./(mi+mi);

lambda_De = sqrt(epsilon0_const.*Te./e_const.^2./ne);
vthe = sqrt(2.*Te./me);
vthi = sqrt(2.*Ti./mi);

lnAei = log(lambda_De./(e_const.^2./4./pi./epsilon0_const./mrei./vthe.^2));
lnAee = log(lambda_De./(e_const.^2./4./pi./epsilon0_const./mree./vthe.^2));
lnAii = log(lambda_De./(e_const.^2./4./pi./epsilon0_const./mrii./vthi.^2));

vei = 2./3./sqrt(2.*pi).*ne.*Z.^2.*e_const.^4./(4.*pi.*epsilon0_const).^2.*4.*pi./sqrt(me.*Te.^3).*lnAei;
vee = 1./3./sqrt(pi).*ne.*e_const.^4./(4.*pi.*epsilon0_const).^2.*4.*pi./sqrt(me.*Te.^3).*lnAee;
vie = me./mi.*vei;
vii = 1./3./sqrt(pi).*ne.*e_const.^4./(4.*pi.*epsilon0_const).^2.*4.*pi./sqrt(mi.*Ti.^3).*lnAii;
end


