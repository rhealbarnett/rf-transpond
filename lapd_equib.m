%------------------------------------------------------------------%
% call equilibrium profiles for the lapd transport solve
% rlbarnett c3149416 190304 
%------------------------------------------------------------------%


% equib = load('../../lapd_numdata/matlab/equibhe_8m_refined.mat');
% equib = load('/Volumes/DATA/LAPD/matlab/lapd_equib_refined.mat');
% equib = load('/Volumes/DATA/LAPD/matlab/lapd_equib_superrefined.mat');
equib = load('/Users/rhealbarnett/Downloads/lapd_equib_superrefined.mat');
% equib = load('C:\Users\c3149416\Documents\lapd_equib_superrefined.mat');
% equib = load('C:\Users\c3149416\Documents\lapd_equib_refined.mat');

vxax = equib.vxax;
nxax = equib.nxax;
vdx = equib.vdx;
ndx = equib.ndx;
dt = equib.dt;
Te = equib.Te;
Ti = equib.Ti;
npts = equib.npts;
tmax = equib.tmax;
n_source = equib.n_source;
vx_new = equib.vx_new;
n_new = equib.n_new; 
xmin = equib.xmin;
xmax = equib.xmax;
nu = equib.nu;

%%

plots = 0;

if plots
    figure(1)
    plot(nxax,n_new,'.-b');
    hold on
else
end

%%

Nmax = 1.0e17;
fact = Nmax/max(n_new);
n_new = n_new*fact;
n_source = n_source*fact;
n_init = n_new;

%%

if plots 
    figure(1); hold on
    plot(nxax,n_new);
    xlabel('Position (m)')
    ylabel('Density (m^{-3})')
    xlim([xmin xmax])
    hold on
else
end

%%

refine = 0;

if refine

    n_source(1,1) = interp1([nxax(3), nxax(2)], [n_source(3), n_source(2)],...
            nxax(1),'linear','extrap'); 
    n_source(1,end) = interp1([nxax(npts-2), nxax(npts-1)], [n_source(npts-2), n_source(npts-1)],...
            nxax(npts),'linear','extrap');

    figure(20);
    plot(nxax(2:npts-1),n_source(2:npts-1))
    hold on

    NP = 2048;
    variable = 1;

    if variable

        % location of sign change
        xc = (xmax - 0)/2.0;

        %'strength' of grid refinement.
        % sign also indicates whether refinement is in the centre or at the
        % boundaries
        A = -4.5;

        % set up the unit spaced parameter, xi, that the grid is a function of
        smax = 1.0;
        smin = 0.0;
        s = linspace(smin,smax,NP-1); 

        % calculate the x values from xi
        x = xc*(1.0 - tanh(A*(1.0 - 2.0*s))./tanh(A));

        Vxax = x - xmax;
        Vdx = (Vxax(2:end) - Vxax(1:end-1));
        
        Vdx = [Vdx, fliplr(Vdx)];
        Vxax = [Vxax, fliplr(-1*Vxax(1:end-1))];

        NP = length(Vxax) + 1;
        Nxax = zeros(1,NP);

        Nxax(1) = Vxax(1) - 0.5*Vdx(1);
        Nxax(end) = Vxax(end) + 0.5*Vdx(end);
        Nxax(2:end-1) = (Vxax(2:end) + Vxax(1:end-1))/2.0;
        Ndx = Nxax(2:end) - Nxax(1:end-1);

    end

    n_new = interp1(nxax,n_new,Nxax,'linear');
    vx_new = interp1(vxax,vx_new,Vxax,'linear');
    n_source = interp1(nxax,n_source,Nxax,'linear');

    n_source(1,1) = 0.0;
    n_source(1,end) = 0.0;

    vxax = Vxax;
    nxax = Nxax;
    ndx = Ndx;
    vdx = Vdx;
    npts = NP;

    if plots
    
        figure(20)
        plot(nxax(2:npts-1),n_source(2:npts-1))
        hold off

        figure(1);
        hold on
        plot(nxax, n_new)
        hold off
        
    end
    
end


%%

lapd_params;
m = real(mhe(1));

cs = sqrt((Te+Ti)*abs(e)/m);
LuBC = -cs;
RuBC = cs;
 
dt = 0.99*min(ndx)/cs;

% source_mult = 37000;
period = 1.0/freq;
% tmax = 100*period;
% tmax = 5.0e-3;
save_time = period/10.0;
% nmax = round(tmax/dt);
nmax = 10;
save_iter = round(save_time/dt);

n_new_uni = interp1(nxax,n_new,xax,'linear');

[om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new_uni,m_s,om,eps0,npts,1);
[A,rf_e,rf_ex,rf_ey,rf_ez] = wave_sol(xax,ky,kx,k0,...
    om,mu0,cpdt,source*1.0e-3,0,1,0);

% rf_ez = zeros(1,npts);
Efield = abs(rf_ez).^2;
Emag = max(abs(sqrt(Efield)));


%%

if plots
    figure(10) 
    plot(xax,sqrt(Efield),'-o')
    hold on
else
end

Efield = interp1(xax,Efield,vxax,'linear');

if plots
    figure(10)
    plot(vxax,sqrt(Efield),'*-')
    hold off
else
end

%%

Ex = interp1(xax,rf_ex,vxax,'linear');
Ey = interp1(xax,rf_ey,vxax,'linear');
Ez = interp1(xax,rf_ez,vxax,'linear');
pond = pond_source({'total',0},{Ex,Ey,Ez},m_s,q_s,om_c,om,vdx,1,{1,vxax});

vx_mat = sparse(nmax,npts-1);
n_mat = sparse(nmax,npts);
pressure_mat = sparse(nmax,npts-2);

vx_mat(1,:) = vx_new;
n_mat(1,:) = n_new;

nA = sparse(npts,npts);
nI = sparse(eye(npts,npts));
vx_pos = sparse(npts-1,npts-1);
vx_neg = sparse(npts-1,npts-1);
vx_diff = sparse(npts-1,npts-1);
vx_I = sparse(eye(npts-1,npts-1));