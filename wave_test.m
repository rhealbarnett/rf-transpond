%-----------------------------------------%
% Run wave solver for lapd-like 
% parameters. 
%-----------------------------------------%

kirallie_waveinput;
[om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,npts,{0,damp_len});
% kz_dispersion_left = dispersion(1,s_arr(1),d_arr(1),p_arr(1),...
%         om,n_refrac,n_new(1),1,0);
% kz_dispersion_right = dispersion(1,s_arr(npts),d_arr(npts),p_arr(npts),...
%         om,n_refrac,n_new(npts),1,0);
[A,rf_e,rf_ex,rf_ey,rf_ez] = wave_sol(zax,ky,kx,k0,...
        om,mu0,cpdt,source,0,0,0,0);
% wave_sol_plots
poyn = poynting(rf_ex, rf_ey, rf_ez, kx, ky, zax, om);
