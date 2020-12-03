%-----------------------------------------%
% Calculate plasma frequency              %
%                                         %
% rlbarnett c3149416, 281018              %
%-----------------------------------------%

function [ans] = plasma_freq(q,n,m,eps0)

    ans = sqrt(n*q^2./(m*eps0));

end
