%-----------------------------------------%
% Run wave solver for lapd-like 
% parameters. 
%-----------------------------------------%

lapd_params;
[om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,npts,1);
[A,rf_e,rf_ex,rf_ey,rf_ez,diss_pow] = wave_sol(xax,ky,kx,k0,...
om,mu0,cpdt,source,0,1,0);
wave_sol_plots