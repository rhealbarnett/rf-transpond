function [ans] = mms_source_mom(om,ux,kux,vxax,dt,ii,nu,u,nxax,knx,nx,n,npts)
    dudt = -om*ux*sin(kux*vxax.^2 + om*dt*ii);
    dudx = -2.0*kux*ux*vxax.*sin(kux*vxax.^2 + om*dt*ii);
    d2udx = -2.0*kux*ux*sin(kux*vxax.^2 + om*dt*ii) -...
        4.0*kux^2*ux*vxax.^2.*cos(kux*vxax.^2 + om*dt*ii);
    dndx = 2.0*knx*nx*nxax.*cos(knx*nxax.^2 + om*dt*ii);
    dndx = interp1(nxax,dndx,vxax);
%     dndx = zeros(1,npts-1);
%     const = constants();
%     e = const.e;
    ans = dudt + u.*dudx - nu*d2udx + ((0.5 + 0.5)./avg(n,npts)).*dndx;
end