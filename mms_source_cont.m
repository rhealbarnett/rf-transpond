function [ans] = mms_source_cont(om,nx,knx,nxax,dt,ii,u,kux,ux,vxax,n)
    dndt = nx*om*cos(knx*nxax.^2 + om*dt*ii);
    dnudx = 2*knx*nx*nxax.*cos(knx*nxax.^2 + om*dt*ii).*u -...
        2*kux*ux*vxax.*sin(kux*vxax.^2 + om*dt*ii).*n;
%     dnudx = -2.0*nx*knx*ux*nxax.*sin(knx*nxax.^2).*cos(kux*vxax.^2) -...
%         2.0*ux*nx*kux*vxax.*cos(knx*nxax.^2)*sin(kux*vxax.^2);
    ans = dndt + dnudx;
end