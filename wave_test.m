%-----------------------------------------%
% Run wave solver for lapd-like 
% parameters. 
%-----------------------------------------%

lapd_params;
[om_c,om_p,cpdt,s_arr,d_arr,p_arr] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,npts);
[ey_source,R] = e_source(xax,scale_fact,source_dist);
[A,source,rf_e,rf_ex,rf_ey,rf_ez] = wave_sol(xax,ky,kx,k0,...
    om,mu0,cpdt,sig,ey_source,R,0,1);
wave_sol_plots