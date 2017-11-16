%------------------------------------------%
% calculate static electric field          %
% E = -grad*static_pot                     %
% rlbarnett c3149416 071117                %
%------------------------------------------%

static_ey = -lamby*static_pot;
static_ez = -lambz*static_pot;

% static_ex = zeros(npts,1);

% for ii=2:(npts-1)
%     static_ex(ii) = (static_pot(ii-1) - static_pot(ii+1))/(2.0*dx);
% end

% static_ex(1) = (-3.0*static_pot(1) + 4.0*static_pot(2) - static_pot(3))/(2.0*dx);
% static_ex(npts) = (3.0*static_pot(npts) - 4.0*static_pot(npts-1) + static_pot(npts-2))/(2.0*dx);

static_ex = gradient(static_pot,dx);