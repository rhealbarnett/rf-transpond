function tests = rf_transp_tests

    tests = functiontests(localfunctions);

end

function [] = mms_test(testCase)

    [linf_arr,ltwo_arr,ratio_inf_arr,ratio_two_arr,...
    oo_inf_arr,oo_two_arr] = run_mms(0,1,1,0);

end