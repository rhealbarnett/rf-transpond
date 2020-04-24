function d_pow = diss_pow(rf_ex,rf_ey,rf_ez,sig)

    exj_dot = zeros(1,npts);
    eyj_dot = zeros(1,npts);
    ezj_dot = zeros(1,npts);
    
    for ii=1:npts
        
        exj_dot(1,ii) = rf_ex(1,ii)*(conj(sig(1,1,ii)*conj(rf_ex(1,ii)) +...
            conj(sig(1,2,ii))*conj(rf_ey(1,ii)) + conj(sig(1,3,ii))*conj(rf_ez(1,ii)))) +...
            conj(rf_ex(1,ii))*(sig(1,1,ii)*rf_ex(1,ii) + sig(1,2,ii)*rf_ey(1,ii) +...
            sig(1,3,ii)*rf_ez(1,ii));
        eyj_dot(1,ii) = rf_ey(1,ii)*(conj(sig(2,1,ii)*conj(rf_ex(1,ii)) +...
            conj(sig(2,2,ii))*conj(rf_ey(1,ii)) + conj(sig(2,3,ii))*conj(rf_ez(1,ii)))) +...
            conj(rf_ey(1,ii))*(sig(2,1,ii)*rf_ex(1,ii) + sig(2,2,ii)*rf_ey(1,ii) +...
            sig(2,3,ii)*rf_ez(1,ii));
        ezj_dot(1,ii) = rf_ez(1,ii)*(conj(sig(3,1,ii)*conj(rf_ex(1,ii)) +...
            conj(sig(3,2,ii))*conj(rf_ey(1,ii)) + conj(sig(3,3,ii))*conj(rf_ez(1,ii)))) +...
            conj(rf_ez(1,ii))*(sig(3,1,ii)*rf_ex(1,ii) + sig(3,2,ii)*rf_ey(1,ii) +...
            sig(3,3,ii)*rf_ez(1,ii));
        
    end
    
    d_powx = (1.0/4.0)*real(trapz(ax,exj_dot));
    d_powy = (1.0/4.0)*real(trapz(ax,eyj_dot));
    d_powz = (1.0/4.0)*real(trapz(ax,ezj_dot));

    d_pow = d_powx + d_powy + d_powz;
    
end