function [ans] = grad2(n,dx,npts)
    cen_diff = (n(3:npts) - n(1:npts-2))/(2.0*dx);
    fwd_diff = (-3*n(1) + 4*n(2) - n(3))/(2.0*dx);
    bwd_diff = (3*n(npts) - 4*n(npts-1) + n(npts-2))/(2.0*dx);
    ans = [fwd_diff, cen_diff, bwd_diff];
end