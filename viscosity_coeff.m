%--------------------------------------------------------------------------------
% Calculate viscosity coefficient        
% eta_i = 0.96*n_i*tau_i*T_i             
% From collisional transport in magnetised plasmas, helander & sigmar  
% rlbarnett c3149416 021018              
%--------------------------------------------------------------------------------
% Having trouble calculating the viscosity coefficient in the parallel
% momentum equation. This book mentions that all calcs are in SI units
% apart from temperature, which is in joules. 
% I think I'm missing something... Because the viscosity comes out to be
% way too small. I did a unit check and it seems like the kinematic
% viscosity is the term we want, given the units.
%--------------------------------------------------------------------------------

%%
% constants;
% 
% mi = const.mi;
% e = const.e;
% eps0 = const.eps0;

mi = 1.67e-27;
e = 1.6022e-19;
eps0 = 8.85e-12;

% temp 5 eV
Ti = 5.0;
% density m^-3
ni = 1.0e18;
% hydrogen
Z = 1.0;
% coulomb logarithm is generally between 10 and 20
co_lo = 15;


% tau_ii taken from page 5 of helander & sigmar
tau_ii = ((12.0*(pi^(3/2)))/(sqrt(2.0)))*...
    ((sqrt(mi)*(Ti*e)^(3/2)*eps0^2)/(ni*Z^4*e^4*co_lo));
% mentions that it is common to multiply this by sqrt(2), per Braginskii
tau_i = sqrt(2)*tau_ii;

% dynamic viscosity kg*(ms)^-1
eta_i = 0.96*ni*(Ti*e)*tau_i;

% kinematic viscosity (has correct units for the momentum equation I am
% using kg*(m^2*s^-1)
% BUT seems waaaaaaaay too small?? 
visc_i = eta_i/ni;

