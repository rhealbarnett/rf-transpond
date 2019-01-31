% npts = 4096;
pow = 3;

for kk=1:3
    
    npts = 2^pow;
    
    fprintf('Iteration %d\n',kk)
    fprintf('Number of points %d\n',npts)
    
    transport_1d
    
    ltwon_arr(1,kk) = l_twon;
    linfn_arr(1,kk) = l_infn;
    ltwou_arr(1,kk) = l_twou;
    linfu_arr(1,kk) = l_infu;
    
    pow = pow + 1;
    
end