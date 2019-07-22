%%
%-------------------------------------------------------------------------%
% CALCULATE DENSITY SOURCE                                                %
%-------------------------------------------------------------------------%
% 
% rate coefficient (constant)
rate_coeff = (1.0e-14);
% approx size of non-zero portion of neutral profile (1/4 domain)
decay_loc = xmax - 0.0001*xmax;
% decay_loc = xmax - 0.2*xmax;
a = find(nxax >= decay_loc);
decay_index = npts - a(1);

% decay_index = round((npts)/4);
% calculate shape of neutral profile
cosax = linspace(pi,2*pi,decay_index);
% max neutral value (at wall)
neut_max = (1.0e18);

% initialise and fill neutral density array
n_neut = zeros(1,npts);
n_neut(end-decay_index+1:end) = neut_max*((cos(cosax) + 1)/2);
% n_neut(1:end-decay_index+1) = n_neut(end-decay_index+2);
n_neut(end/2 + 1:end-decay_index+1) = neut_max*((cos(pi) + 1)/2);
n_neut(1,1:end/2) = fliplr(n_neut(end/2 + 1:end));
n_neut = interp1(linspace(xmin-0.5*dx,xmax+0.5*dx,npts),n_neut,nxax,'linear');

% calculate density source
n_source = (n_new.*(n_neut)*(rate_coeff));

% interpolate source onto velocity grid
source_avg = interp1(nxax,n_source,vxax);
% calculate integral of density source over the velocity grid
source_int = trapz(vxax,source_avg);
% source_int = trapz(nxax,n_source);
% calculate flux at the rh boundary (wall)
% rflux = n_avg(end)*vx_new(end);
rflux = n_avg(end)*vx_new(end);
lflux = n_avg(1)*vx_new(1);
% calculate the constant multiplier to match density out = in
ns_mult = (rflux-lflux)/source_int;
% multiply n0(x)n(x,t) by the constant calculated in previous step
n_source = (n_source*ns_mult)*1.0e-2;
nv_source = source_avg*ns_mult;
n_source(1,1) = 0.0; n_source(1,end) = 0.0;