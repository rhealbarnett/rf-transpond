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

    fprintf(['=======================================',...
        '\n Running time dependent transport MMS\n',...
        '=======================================\n'])

    [linf_arr,ltwo_arr,ratio_inf_arr,ratio_two_arr,...
    oo_inf_arr,oo_two_arr,arr] = run_mms(0,1,1,0);

    p = polyfit(log10(arr),log10(ltwo_arr(1,:)),1);

    expSlope = 1.0;
    actSlope = p(1);
    
    verifyEqual(testCase,actSlope,expSlope,'RelTol',1e-1)

end