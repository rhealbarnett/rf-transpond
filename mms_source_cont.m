function [ans] = mms_source_cont(om,nx,knx,nxax,dt,ii,u,kux,ux,vxax,n)
    dndt = nx*om*cos(knx*nxax.^2 + om*dt*ii);
    dnudx = 2*knx*nx*nxax.*cos(knx*nxax.^2 + om*dt*ii).*u -...
        2*kux*ux.*vxax.*sin(kux*vxax.^2 + om*dt*ii).*n;
    ans = dndt + dnudx;
end