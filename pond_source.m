%------------------------------------------------------------------%
% Function calculates ponderomotive acceleration
% rlbarnett c3149416 updated 191127.
%
% INPUTS
% m: mass of the particle being accelerated.
% mp: proton mass from the velocity update equation.
% omega: driving frequency for the fast E field.
% q: charge of the particle being accelerated.
% Efield: absolute value of the fast electric field squared, the quantity
%       to be differentiated.
% dx: array of dx values. 
% varargin: option to include cyclotron frequency required for calculation
%       of the ion ponderomotive acceleration.
%
% OUPUTS
% ans: the ponderomotive acceleration for input into an eq of motion, the 
%       velocity equation.
%
%------------------------------------------------------------------%

function [ans] = pond_source(m,mp,omega,q,Efield,dx,varargin)
    for ii = 2:length(dx)
        Ediff(1,ii-1) = -(dx(1,ii)/(dx(1,ii-1)*(dx(1,ii) + dx(1,ii-1))))*Efield(1,ii-1) +...
            ((dx(1,ii) - dx(1,ii-1))/(dx(1,ii)*dx(1,ii-1)))*Efield(1,ii) +...
            (dx(1,ii-1)/(dx(1,ii)*(dx(1,ii) + dx(1,ii-1))))*Efield(1,ii+1);
    end
    
    if numel(varargin)==0
        pond_const = (1.0/4.0)*((q^2)/(m*omega^2));
    else
        pond_const = (1.0/4.0)*((q^2)/(m*(omega^2 - varargin{1}^2)));
    end
    
    ans = (1.0/mp)*pond_const*Ediff;
end