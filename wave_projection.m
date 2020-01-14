%-------------------------------------------------------------------------%
% Projection of wave solution into ignorable directions. 
% rlbarnett c3149416 200109
%-------------------------------------------------------------------------%

function [Ex_x, Ey_x, Ez_x, Ex_y, Ey_y, Ez_y] = wave_projection(xax,yax,kx,ky,rf_ex,rf_ey,rf_ez,npts)

    nptsx = numel(xax);
    nptsy = numel(yax);

    Ex_x = zeros(npts,nptsx);
    Ex_y = zeros(npts,nptsy);
    Ey_x = zeros(npts,nptsx);
    Ey_y = zeros(npts,nptsy);
    Ez_x = zeros(npts,nptsx);
    Ez_y = zeros(npts,nptsy);
    
    for jj=1:npts
        for ii=1:nptsx

            if xax(1,ii)<0
                xax(1,ii) = -1*xax(1,ii);
            end

            Ex_x(jj,ii) = real(rf_ex(1,jj)*(exp(1i*kx(1,ii)*xax(1,ii))));
            Ey_x(jj,ii) = real(rf_ey(1,jj)*(exp(1i*kx(1,ii)*xax(1,ii))));
            Ez_x(jj,ii) = real(rf_ez(1,jj)*(exp(1i*kx(1,ii)*xax(1,ii))));

        end

        for ii=1:nptsy

            if yax(1,ii)<0
                yax(1,ii) = -1*yax(1,ii);
            end
            
            Ex_y(jj,ii) = real(rf_ex(1,jj)*(exp(1i*ky(1,ii)*yax(1,ii))));
            Ey_y(jj,ii) = real(rf_ey(1,jj)*(exp(1i*ky(1,ii)*yax(1,ii))));
            Ez_y(jj,ii) = real(rf_ez(1,jj)*(exp(1i*ky(1,ii)*yax(1,ii))));

        end
    end

end
