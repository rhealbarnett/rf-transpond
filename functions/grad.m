function [ans] = grad(n,dx,npts)
    ans = (n(2:npts) - n(1:npts-1))./dx;
end