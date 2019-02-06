% npts = 4096;
pow = 4;

for kk=1:6
    
    npts = 2^pow;
    
    fprintf('Iteration %d\n', kk)
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
    
end