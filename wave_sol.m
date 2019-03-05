%------------------------------------------------------------------%
% time independent wave solver                                   
% V X V X E = k0^2.K.E           
% rlbarnett c3149416 210917      
%------------------------------------------------------------------%
% edited rlbarnett 291018
% changed to a function
% need to run dielec_tens function to get cpdt before calling wave_sol
%------------------------------------------------------------------%


function [A,source,rf_e,rf_ex,rf_ey,rf_ez,diss_pow] = wave_sol(xax,ky,kz,k0,...
    om,mu0,cpdt,source_width,source_loc,MMS)

    npts = length(xax);
    xmax = xax(1,end);
    xmin = xax(1,1);
    dx = (xmax - xmin)/(npts-1);
    A = sparse(3*npts, 3*npts);
    ii = 4;
    kk = 2;

    for eq1=4:3:3*(npts-1)
        
        eq2 = eq1 + 1;
        eq3 = eq2 + 1;

        iiex = ii;
        iiey = ii+1;
        iiez = ii+2;
        iiexm = iiex - 3;
        iieym = iiey - 3;
        iiezm = iiez - 3;
        iiexp = iiex + 3;
        iieyp = iiey + 3;
        iiezp = iiez + 3;        

        %--
        % fill matrix
        A(eq1,iiexm) = 0.0;
        A(eq1,iieym) = -1i*ky/(2.0*dx);
        A(eq1,iiezm) = -1i*kz/(2.0*dx);
        A(eq1,iiex) = (ky^2 + kz^2 - k0^2*cpdt(1,1,kk));
        A(eq1,iiey) = -k0^2*cpdt(1,2,kk);
        A(eq1,iiez) = -k0^2*cpdt(1,3,kk);
        A(eq1,iiexp) = 0.0;
        A(eq1,iieyp) = 1i*ky/(2.0*dx);
        A(eq1,iiezp) = 1i*kz/(2.0*dx);

        A(eq2,iiexm) = -1i*ky/(2.0*dx);
        A(eq2,iieym) = -1.0/(dx^2);
        A(eq2,iiezm) = 0.0;
        A(eq2,iiex) = -k0^2*cpdt(2,1,kk);
        A(eq2,iiey) = (kz^2 - k0^2*cpdt(2,2,kk)) + 2.0/(dx^2);
        A(eq2,iiez) = -ky*kz - k0^2*cpdt(2,3,kk);
        A(eq2,iiexp) = 1i*ky/(2.0*dx);
        A(eq2,iieyp) = -1.0/(dx^2);
        A(eq2,iiezp) = 0.0;

        A(eq3,iiexm) = -1i*kz/(2.0*dx);
        A(eq3,iieym) = 0.0;
        A(eq3,iiezm) = -1.0/(dx^2);
        A(eq3,iiex) = -k0^2*cpdt(3,1,kk);
        A(eq3,iiey) = -ky*kz - k0^2*cpdt(3,2,kk);
        A(eq3,iiez) = (ky^2 - k0^2*cpdt(3,3,kk)) + 2.0/(dx^2);
        A(eq3,iiexp) = 1i*kz/(2.0*dx);
        A(eq3,iieyp) = 0.0;
        A(eq3,iiezp) = -1.0/(dx^2);

        ii = ii + 3;
        kk = kk + 1;

    end

    %--
    % metallic wall BC
    A(1,1) = 1.0;
    A(2,2) = 1.0;
    A(3,3) = 1.0;
    A(end-2,end-2) = 1.0;
    A(end-1,end-1) = 1.0;
    A(end,end) = 1.0;
    
    %--
    % set up rhs vector (current source term)
    rhs = zeros(3*npts,1);
    mult = 1.0/sqrt(2.0*pi*source_width);
    source = mult*exp(-(xax - source_loc).^2/(2.0*source_width^2));
    source = source / max(source);
    source = source*10000;
    rhs(1:3:3*npts,1) = 1i*om*mu0*source';
    rhs(2:3:3*npts,1) = 1i*om*mu0*source';
    rhs(3:3:3*npts,1) = 1i*om*mu0*source';
    
    if MMS
            
            Ex = 1.0;
            Ey = 1.0;
            Ez = 0.01;
           
            kx = 40.;
            
            ex_solx = Ex*exp(1i*kx*xax);
            ex_soly = Ey*exp(1i*kx*xax);
            ex_solz = Ez*exp(1i*kx*xax);
            
            ex_sol = zeros(1,3*npts);
            ex_sol(1:3:3*npts) = ex_solx;
            ex_sol(2:3:3*npts) = ex_soly;
            ex_sol(3:3:3*npts) = ex_solz;
            
            source_x = zeros(1,npts);
            source_y = zeros(1,npts);
            source_z = zeros(1,npts);
            
        for jj=1:npts
            
            source_x(1,jj) = 1i*ky*(1i*kx*ex_soly(1,jj)) +...
                (ky^2 + kz^2)*ex_solx(1,jj) + 1i*kz*(1i*kx*ex_solz(1,jj)) -...
                k0^2*(cpdt(1,1,jj).*ex_solx(1,jj) + cpdt(1,2,jj).*ex_soly(1,jj) +...
                cpdt(1,3,jj).*ex_solz(1,jj)) + 1i*om*mu0*source(1,jj);
            source_y(1,jj) = 1i*ky*(1i*kx*ex_solx(1,jj)) +...
                (kx^2 - kz^2)*ex_soly(1,jj) - ky*kz*ex_solz(1,jj) -...
                k0^2*(cpdt(2,1,jj).*ex_solx(1,jj) + cpdt(2,2,jj).*ex_soly(1,jj) +...
                cpdt(2,3,jj).*ex_solz(1,jj)) + 1i*om*mu0*source(1,jj);
            source_z(1,jj) = 1i*kz*(1i*kx*ex_solx(1,jj)) +...
                (kx^2 + ky^2)*ex_solz(1,jj) - ky*kz*ex_soly(1,jj) -...
                k0^2*(cpdt(3,1,jj).*ex_solx(1,jj) + cpdt(3,2,jj).*ex_soly(1,jj) +...
                cpdt(3,3,jj).*ex_solz(1,jj)) + 1i*om*mu0*source(1,jj);
            
        end
        
        source_mms = zeros(1,3*npts);
        source_mms(1:3:3*npts) = source_x;
        source_mms(2:3:3*npts) = source_y;
        source_mms(3:3:3*npts) = source_z;
        
        source_mms(1,1) = Ex*exp(1i*kx*xmin); 
        source_mms(1,2) = Ey*exp(1i*kx*xmin);
        source_mms(1,3) = Ez*exp(1i*kx*xmin);
        source_mms(1,end-2) = Ex*exp(1i*kx*xmax); 
        source_mms(1,end-1) = Ey*exp(1i*kx*xmax);
        source_mms(1,end) = Ez*exp(1i*kx*xmax);
        
        rf_e = A\(source_mms' + rhs);
        
        rf_e = rf_e';
        
    else

        rf_e = (A)\rhs;
        rf_e = rf_e';
        
    end

    rf_ex = rf_e(1,1:3:3*npts);
    rf_ey = rf_e(1,2:3:3*npts);
    rf_ez = rf_e(1,3:3:3*npts);
    
    exj_dot = rf_ex.*conj(source);
    eyj_dot = rf_ey.*conj(source);
    ezj_dot = rf_ez.*conj(source);
    
    diss_powx = (1.0/2.0)*real(trapz(xax,exj_dot));
    diss_powy = (1.0/2.0)*real(trapz(xax,eyj_dot));
    diss_powz = (1.0/2.0)*real(trapz(xax,ezj_dot));
    
    diss_pow = [diss_powx, diss_powy, diss_powz];

end

