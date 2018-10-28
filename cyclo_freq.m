%-----------------------------------------%
% Calculate cyclotron frequency           %
%                                         %
% rlbarnett c3149416, 281018              %
%-----------------------------------------%

function [ans] = cyclo_freq(q,B0,m)

    ans = q*B0./m;

end