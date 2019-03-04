pow = 4;

for kk=1:9
    
    npts = 2^pow;
    
    fprintf('iteration %d\n', kk)
    fprintf('Number of points %d\n', npts)
    
    lapd_params
    [om_c,om_p,cpdt,s_arr,d_arr,p_arr] = dielec_tens(charge,B0,n_new,m,om,eps0,npts);
    [A,source,rf_e,rf_ex,rf_ey,rf_ez,ex_sol] = wave_sol(xax,ky,kz,k0,om,mu0,cpdt,...
    source_width,xmax/2,1);

    l_inf = norm(ex_sol - rf_e, Inf);
    l_two = rms(ex_sol - rf_e);

    ltwo_arr(1,kk) = l_two;
    linf_arr(1,kk) = l_inf;
    
    pow = pow + 1;
    npts_arr(1,kk) = npts;
    
end

dx_arr = 1.0./npts_arr;

ratio_inf = linf_arr(1:kk-1)./linf_arr(2:kk);
ratio_two = ltwo_arr(1:kk-1)./ltwo_arr(2:kk);
oo_inf = log(ratio_inf)/log(2);
oo_two = log(ratio_two)/log(2);