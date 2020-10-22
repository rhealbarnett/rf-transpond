function [ans] = pressure_source_col(n,q,T,m,npts,dx)
    ans = -((T)*q./(m*n)).*(grad2(n,dx,npts));
end