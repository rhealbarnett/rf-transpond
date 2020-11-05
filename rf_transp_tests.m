function tests = rf_transp_tests

    tests = functiontests(localfunctions);

end

function [] = mms_ss_test(testCase)

    fprintf(['=======================================',...
        '\n Running steady state transport MMS\n',...
        '=======================================\n'])

    [linf_arr,ltwo_arr,ratio_inf_arr,ratio_two_arr,...
    oo_inf_arr,oo_two_arr,arr] = run_mms(1,0,1,0);

    p = polyfit(log10(arr),log10(ltwo_arr(1,:)),1);

    expSlope = -1.0;
    actSlope = p(1);
    
    verifyEqual(testCase,actSlope,expSlope,'RelTol',1e-1)

end


function [] = mms_td_test(testCase)

    fprintf(['\n=======================================',...
        '\n Running time dependent transport MMS\n',...
        '=======================================\n'])

    [linf_arr,ltwo_arr,ratio_inf_arr,ratio_two_arr,...
    oo_inf_arr,oo_two_arr,arr] = run_mms(0,1,1,0);

    p = polyfit(log10(arr),log10(ltwo_arr(1,:)),1);

    expSlope = 1.0;
    actSlope = p(1);
    
    verifyEqual(testCase,actSlope,expSlope,'RelTol',1e-1)

end

function[] = dispersion_test(testCase)

    wave_verification;
    
    [om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,npts,{0,''});
    kz_dispersion = dispersion(npts,s_arr,d_arr,p_arr,om,n_refrac,n_new,1,0,ky);
    
    expkz = real(kz_dispersion(1,:));
   
    kz_spec_density = zeros(npts,npts);

    count = 1;

    for ii=1:npts

        density = n_new(1,ii)*ones(1,npts);

        [om_c,om_p,cpdt,s_arr,d_arr,p_arr,sig] = dielec_tens(q_s,B0,density,m_s,om,eps0,npts,{1,damp_len});
        [A,rf_e,rf_ex,rf_ey,rf_ez] = wave_sol(zax,ky,kx,k0,...
        om,mu0,cpdt,source,0,1,1,0);

        [kz_spec, k_ax, phase, dk] = fft_kz(dz,npts,rf_ex,rf_ey,rf_ez,0);

        kz_spec_density(count,:) = kz_spec(:,3);
        
        ind_kz = find(kz_spec(1:npts/2,2)==max(kz_spec(1:npts/2,2)));
        actkz(1,count) = k_ax(ind_kz);

        count = count + 1;

    end
    
    verifyEqual(testCase,actkz,expkz,'RelTol',dk)

end