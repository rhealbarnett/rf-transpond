function [ans] = pressure_source_stag(n,q,Te,Ti,m,npts,dx)
    ans = -((Te + Ti)*q./(m*avg(n,npts))).*(grad(n,dx,npts));
end