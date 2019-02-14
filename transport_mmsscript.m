% npts = 4096;
pow = 4;

for kk=1:9
    
    npts = 2^pow;
    
    fprintf('kkation %d\n', kk)
    fprintf('Number of points %d\n', npts)
    
    transport_1d
    
    fprintf('Minimum dx %d\n', min(ndx))
    fprintf('Maximum dx %d\n', max(ndx))
    
%     ltwon_arr{kk} = l_twon;
%     linfn_arr{kk} = l_infn;
%     ltwou_arr{kk} = l_twou;
%     linfu_arr{kk} = l_infu;
    ltwon_arr(1,kk) = l_twon;
    linfn_arr(1,kk) = l_infn;
    ltwou_arr(1,kk) = l_twou;
    linfu_arr(1,kk) = l_infu;
    
    pow = pow + 1;
    npts_arr(1,kk) = npts;
    
end

ratio_inf = linfu_arr(1:kk-1)./linfu_arr(2:kk);
ratio_two = ltwou_arr(1:kk-1)./ltwou_arr(2:kk);
oo_inf = log(ratio_inf)/log(2);
oo_two = log(ratio_two)/log(2);