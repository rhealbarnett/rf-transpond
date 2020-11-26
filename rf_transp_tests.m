function tests = rf_transp_tests

    tests = functiontests(localfunctions);

end

function [] = mms_ss_test(testCase)

    fprintf(['=======================================',...
        '\n Running steady state transport MMS\n',...
        '=======================================\n'])

    [~,ltwo_arr,~,~,~,~,arr] = run_mms(1,0,1,0);

    p = polyfit(log10(arr),log10(ltwo_arr(1,:)),1);

    expSlope = -1.0;
    actSlope = p(1);
    
    verifyEqual(testCase,actSlope,expSlope,'RelTol',1e-1)

end


function [] = mms_td_test(testCase)

    fprintf(['\n=======================================',...
        '\n Running time dependent transport MMS\n',...
        '=======================================\n'])

    [~,ltwo_arr,~,~,~,~,arr] = run_mms(0,1,1,0);

    p = polyfit(log10(arr),log10(ltwo_arr(1,:)),1);

    expSlope = 1.0;
    actSlope = p(1);
    
    verifyEqual(testCase,actSlope,expSlope,'RelTol',1e-1)

end

function[] = dispersion_test(testCase)

    fprintf(['\n=======================================',...
        '\n Running wave solver dispersion verification\n',...
        '=======================================\n'])

    wave_verification;
    
    [~,~,~,s_arr,d_arr,p_arr,~] = dielec_tens(q_s,B0,n_new,m_s,om,eps0,...
        npts,{0,damp_len,dampFac});
    kz_dispersion = dispersion(npts,s_arr,d_arr,p_arr,om,n_refrac,n_new,1,0,ky);
    
    expkz = real(kz_dispersion(1,:));
   
    [actkz, dk] = kz_spectrum(n_new,q_s,m_s,om,npts,damp_len,zax,ky,kx,k0,B0,...
                                source,expkz,1);
    
    verifyEqual(testCase,actkz,expkz,'RelTol',dk)

end