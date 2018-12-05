%-------------------------------------------------------%
% parameter file for coupled transport equations        %
% rlbarnett c3149416 051218                             %
%-------------------------------------------------------%

function params = transport_params()

const = constants();

params.Te = 5.0;
params.Ti = 10.0;
params.nu = 1.0;
params.cs = sqrt((params.Te + params.Ti)*const.e/const.mp);

end