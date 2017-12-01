%------------------------------------------%
% calculate static electric field          %
% E = -grad*static_pot                     %
% rlbarnett c3149416 071117                %
%------------------------------------------%

static_ey = -lamby*static_pot;
static_ez = -lambz*static_pot;

static_ex = -gradient(static_pot,dx);