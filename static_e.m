%------------------------------------------%
% calculate static electric field          %
% E = -grad*static_pot                     %
% rlbarnett c3149416 071117                %
%------------------------------------------%

static_ey = -lamby*static_poty;
static_ez = -lambz*static_potz;

static_ex = zeros(npts,1);

for ii=2,(npts-1);
    static_ex(ii) = (static_potx(ii-1) - static_potx(ii+1))/(2.0*dx);
end

static_ex(1) = (-3.0*static_potx(1) + 4.0*static_potx(2) - static_pot(3))/(2.0*dx);
static_ex(npts) = (3.0*static_potx(npts) - 4.0*static_potx(npts-1) + static_potx(npts-2))/(2.0*dx);