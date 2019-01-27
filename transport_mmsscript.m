npts = 4096;

for kk=1:5
    
    fprintf('Iteration %d\n',kk)
    fprintf('Number of points %d\n',npts)
    
    transport_1d
    
    ltwo_arr(1,kk) = l_two;
    linf_arr(1,kk) = l_inf;
    
    npts = npts/2;
    
end