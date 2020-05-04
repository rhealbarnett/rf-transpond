%-----------------------------------------%
% Run wave solver for lapd-like 
% parameters. 
%-----------------------------------------%

lapd_equib;
[om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,npts,1);
[A,rf_e,rf_ex,rf_ey,rf_ez] = wave_sol(zax,ky,kx,k0,...
om,mu0,cpdt,sig,source,0,1,1);
% Ex = interp1(zax,!rf_ex,vxax,'linear');
% Ey = interp1(zax,rf_ey,vxax,'linear');
% Ez = interp1(zax,rf_ez,vxax,'linear');
% [Ediff, pf_source] = pond_source({'para',0},{Ex,Ey,Ez},m_s,q_s,om_c,om,vdx,1,{1,vxax});
wave_sol_plots

%%

xax = linspace(-2,2,npts/2);
yax = 0.0;%linspace(-2,2,npts/4);

[Ex_x, Ey_x, Ez_x, Ex_y, Ey_y, Ez_y, Emag_x, Emag_y] = wave_projection(xax,yax,kx,ky,rf_ex,...
    rf_ey,rf_ez,npts);

