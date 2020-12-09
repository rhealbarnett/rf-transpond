%------------------------------------------------------------------%
% call equilibrium profiles for the lapd transport solve
% rlbarnett c3149416 190304 
% Updated for CPC manuscript 201021
%------------------------------------------------------------------%

%--
% Load transport equilibrium file
load('inputs/equil_transport_input.mat');

%--
% Scale density and density source to desired max density.
Nmax = 1.0e17;
fact = Nmax/max(n_new);
n_new = n_new*fact;
n_source = n_source*fact;
n_init = n_new;

%--
% Call parameter file
lapd_rfparams;

%--
% Set ion mass 
m = real(mhe(1));

%--
% Calculate sound speed for velocity boundary condition
cs = sqrt((Te+Ti)*abs(e)/m);

%--
% Set let and right velocity boundary conditions
LuBC = -cs;
RuBC = cs;

%--
% Set initial velocity
vx_init = vx_new;
 
%--
% Calculate time step based on smallest dx and cs.
dt = 0.99*min(ndx)/cs;

%--
% Set max time based on
% number of periods, number of iterations nmax.
tmax = 2*period;
nmax = round(tmax/dt);
ii = 0;

%--
% Set save frequency
save_time = period/10.0;
save_iter = round(save_time/dt);

%% ---------------------------------------------------------------- %%
% Initial call to wave solver

%--
% Interpolate density onto a uniform grid
n_new_uni = interp1(nxax,n_new,zax,'linear');

%--
% Calculate dielectric tensor
[om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new_uni,m_s,om,eps0,...
    npts,{1,damp_len,dampFac});

%--
% Calculate initial RF wave fields
[A,rf_e,rf_ex,rf_ey,rf_ez] = rf_wave_sol(zax,ky,kx,k0,...
    om,mu0,cpdt,source,0,1,1);

%--
% Calculate initial Poynting flux
poyn = poynting(rf_ex, rf_ey, rf_ez, kx, ky, zax, om);

%--
% Calculate initial ponderomotive force
rampFac = 1.0e-3;
[Ediff, pf] = pond_source({'total',0},{rf_ex,rf_ey,rf_ez},m_s,...
    q_s,om_c,om,dz,1,{1,damp_len,zax});
pf_inter = sum(pf,1);
pf_inter2 = squeeze(sum(pf_inter,2))';
pf_source = interp1(zax,pf_inter2,vxax,'linear');
pf_source(1,1) = 0.0; pf_source(1,end) = 0.0;

%--
% Save struct.

if sfile
    
    transport.dt = dt;
    transport.n_source = n_source;
    transport.vx_new = vx_new;
    transport.n_new = n_new; 
    transport.cs = cs;
    transport.vxax = vxax;
    transport.nxax = nxax;
    transport.vdx = vdx;
    transport.ndx = ndx;
    transport.Te = Te;
    transport.Ti = Ti;
    transport.npts = npts;
    transport.tmax = tmax;
    transport.zmin = xmin;
    transport.zmax = xmax;
    transport.nu = nu;
    transport.freq = freq;
    transport.period = period;
    transport.B0 = B0;
    transport.rf_ez = rf_ez;
    transport.rf_ey = rf_ey;
    transport.rf_ex = rf_ex;
    transport.zax = zax;
    transport.source = source;
    transport.kx = kx;
    transport.ky = ky;
    transport.Ediff = Ediff;
    transport.pond = pf;
    transport.pond_summed = pf_source;
    transport.poyn = poyn;

    [status,git_hash] = system('git rev-parse HEAD');
    s1 = '# Created from matlab git hash ';
    s2 = git_hash;
    header = [s1 s2];

    transport.header = header;

    filename = strcat('outputs/coupled_results/coupled_transport_0.mat');
    save(filename,'-struct','transport');
end
