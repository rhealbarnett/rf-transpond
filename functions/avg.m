function [ans] = avg(n,npts)
    ans = (n(2:npts) + n(1:npts-1))/2.0;
end