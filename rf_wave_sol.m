%------------------------------------------------------------------%
% time independent wave solver, cartesian coordinates.                                   
% V X V X E = k0^2.K.E           
% rlbarnett c3149416 210917      
%------------------------------------------------------------------%

%------------------------------------------------------------------------%
% INPUTS
% ax: spatial axis.
% ky: wavenumber in y, in m^-1.
% k: other wavenumber, either kx or kz depending on solution domain, in
%    m^-1.
% k0: vacuum wavenumber, in m^-1.
% om: driving frequency, in rad s^-1. 
% mu0: vacuum permeability, in SI units.
% cpdt: cold plasma dielectric tensor, 3X3Xnpts, calculated from 
%       dielec_tens.m function (or otherwise if you're fancy).
% source: current source term defined at each location in space, in Am^-2.
% MMS: switch 1 or 0 to run or not run MMS tests for wave solver.
% para: switch 1 or 0 to solve in parallel (z) or perpendicular (x) domain.
% sparsefill: switch 1 or 0 to use or not use matlab's sparsefilling
%             routines to fill the matrix. 
%
% OUTPUTS
% A: electric field coeficient matrix in the wave equation. 
% rf_e(x,y,z): each component of the rf e field solution. 
%------------------------------------------------------------------------%


function [A,rf_e,rf_ex,rf_ey,rf_ez] = rf_wave_sol(ax,ky,k,k0,...
    om,mu0,cpdt,source,MMS,para,sparsefill)

    npts = length(ax);
    axmax = ax(1,end);
    axmin = ax(1,1);
    h = (axmax - axmin)/(npts);
    A = sparse(3*npts, 3*npts);
    ii = 4;
    kk = 2;
    const = constants();
    
    np_sparse = (3*npts - 6)*9 + 6;

    row = zeros(1,np_sparse);
    column = zeros(1,np_sparse);
    v = zeros(1,np_sparse);
    
    row(1,1) = 1;
    row(1,2) = 2;
    row(1,3) = 3;
    row(1,np_sparse-2) = (3*npts)-2;
    row(1,np_sparse-1) = (3*npts)-1;
    row(1,np_sparse) = (3*npts);
    
    column(1,1) = 1;
    column(1,2) = 2;
    column(1,3) = 3;
    column(1,np_sparse-2) = (3*npts)-2;
    column(1,np_sparse-1) = (3*npts)-1;
    column(1,np_sparse) = (3*npts);
    
    v(1,1) = 1;
    v(1,2) = 1;
    v(1,3) = 1;
    v(1,np_sparse-2) = 1;
    v(1,np_sparse-1) = 1;
    v(1,np_sparse) = 1;   
    
    start = 4;
    stop = 12;
    nnex = 7;
    count = 0;
    
    if para && ~sparsefill

        for eq1=4:3:3*(npts-1)

            eq2 = eq1 + 1;
            eq3 = eq2 + 1;       

            iiex = ii;
            iiey = iiex+1;
            iiez = iiex+2;
            iiexm = iiex - 3;
            iieym = iiey - 3;
            iiezm = iiez - 3;
            iiexp = iiex + 3;
            iieyp = iiey + 3;
            iiezp = iiez + 3; 

            %--
            % fill matrix
            A(eq1,iiexm) = -1.0/h^2;
            A(eq1,iieym) = 0.0;
            A(eq1,iiezm) = -1i*k(1,kk)/(2.0*h);
            A(eq1,iiex) = (ky(1,kk)^2 + (2.0/h^2) - k0^2*cpdt(1,1,kk));
            A(eq1,iiey) = -ky(1,kk)*k(1,kk) - k0^2*cpdt(1,2,kk);
            A(eq1,iiez) = -k0^2*cpdt(1,3,kk);
            A(eq1,iiexp) = -1.0/(h^2);
            A(eq1,iieyp) = 0.0;
            A(eq1,iiezp) = 1i*k(1,kk)/(2.0*h);

            A(eq2,iiexm) = 0.0;
            A(eq2,iieym) = -1.0/(h^2);
            A(eq2,iiezm) = -1i*ky(1,kk)/(2.0*h);
            A(eq2,iiex) = -ky(1,kk)*k(1,kk) - k0^2*cpdt(2,1,kk);
            A(eq2,iiey) = (k(1,kk)^2 - k0^2*cpdt(2,2,kk)) + 2.0/(h^2);
            A(eq2,iiez) =  -k0^2*cpdt(2,3,kk);
            A(eq2,iiexp) = 0.0;
            A(eq2,iieyp) = -1.0/(h^2);
            A(eq2,iiezp) = 1i*ky(1,kk)/(2.0*h);

            A(eq3,iiexm) = -1i*k(1,kk)/(2.0*h);
            A(eq3,iieym) = -1i*ky(1,kk)/(2.0*h);
            A(eq3,iiezm) = 0.0;
            A(eq3,iiex) = -k0^2*cpdt(3,1,kk);
            A(eq3,iiey) = -k0^2*cpdt(3,2,kk);
            A(eq3,iiez) = (ky(1,kk)^2 + k(1,kk)^2 - k0^2*cpdt(3,3,kk));
            A(eq3,iiexp) = 1i*k(1,kk)/(2.0*h);
            A(eq3,iieyp) = 1i*ky(1,kk)/(2.0*h);
            A(eq3,iiezp) = 0.0;
            
            ii = ii + 3;
            kk = kk + 1;

        end
        
        
            
    elseif para && sparsefill
        
        row1 = 4;
        row2 = row1 + 1;
        row3 = row2 + 1;
         
        for gg = 1:npts-2
            
            nney = nnex+1;
            nnez = nnex+2;
            nnexm = nnex - 3;
            nneym = nney - 3;
            nnezm = nnez - 3;
            nnexp = nnex + 3;
            nneyp = nney + 3;
            nnezp = nnez + 3; 

            arr = linspace(count*3+1,count*3+9,9);

            row(1,start:stop) = row1;
            row(1,start+9:stop+9) = row2;
            row(1,start+18:stop+18) = row3;
            column(1,start:stop) = arr;
            column(1,start+9:stop+9) = arr;
            column(1,start+18:stop+18) = arr;

            v(1,nnexm) = -1.0/h^2;
            v(1,nneym) = 0.0;
            v(1,nnezm) = -1i*k(1,kk)/(2.0*h);
            v(1,nnex) = (ky(1,kk)^2 + (2.0/h^2) - k0^2*cpdt(1,1,kk));
            v(1,nney) = -ky(1,kk)*k(1,kk) - k0^2*cpdt(1,2,kk);
            v(1,nnez) = -k0^2*cpdt(1,3,kk);
            v(1,nnexp) = -1.0/(h^2);
            v(1,nneyp) = 0.0;
            v(1,nnezp) = 1i*k(1,kk)/(2.0*h);

            v(1,nnexm+9) = 0.0;
            v(1,nneym+9) = -1.0/(h^2);
            v(1,nnezm+9) = -1i*ky(1,kk)/(2.0*h);
            v(1,nnex+9) = -ky(1,kk)*k(1,kk) - k0^2*cpdt(2,1,kk);
            v(1,nney+9) = (k(1,kk)^2 - k0^2*cpdt(2,2,kk)) + 2.0/(h^2);
            v(1,nnez+9) = -k0^2*cpdt(2,3,kk);
            v(1,nnexp+9) = 0.0;
            v(1,nneyp+9) = -1.0/(h^2);
            v(1,nnezp+9) = 1i*ky(1,kk)/(2.0*h);

            v(1,nnexm+18) = -1i*k(1,kk)/(2.0*h);
            v(1,nneym+18) = -1i*ky(1,kk)/(2.0*h);
            v(1,nnezm+18) = 0.0;
            v(1,nnex+18) = -k0^2*cpdt(3,1,kk);
            v(1,nney+18) = -k0^2*cpdt(3,2,kk);
            v(1,nnez+18) = (ky(1,kk)^2 + k(1,kk)^2 - k0^2*cpdt(3,3,kk));
            v(1,nnexp+18) = 1i*k(1,kk)/(2.0*h);
            v(1,nneyp+18) = 1i*ky(1,kk)/(2.0*h);
            v(1,nnezp+18) = 0.0;
        
            start = start + (3*9);
            stop = stop + (3*9);
            nnex = nnex + (3*9);
            count = count + 1;
            row1 = row1 + 3;
            row2 = row1 + 1;
            row3 = row2 + 1;
            kk = kk + 1;
            
        end
        
        
        
    elseif ~para && ~sparsefill



        for eq1=4:3:3*(npts-1)

            eq2 = eq1 + 1;
            eq3 = eq2 + 1;       

            iiex = ii;
            iiey = iiex+1;
            iiez = iiex+2;
            iiexm = iiex - 3;
            iieym = iiey - 3;
            iiezm = iiez - 3;
            iiexp = iiex + 3;
            iieyp = iiey + 3;
            iiezp = iiez + 3;         
        
            %--
            % Perp wave solve
            A(eq1,iiexm) = 0.0;
            A(eq1,iieym) = -1i*ky/(2.0*h);
            A(eq1,iiezm) = -1i*k/(2.0*h);
            A(eq1,iiex) = (ky^2 + k^2 - k0^2*cpdt(1,1,kk));
            A(eq1,iiey) = -k0^2*cpdt(1,2,kk);
            A(eq1,iiez) = -k0^2*cpdt(1,3,kk);
            A(eq1,iiexp) = 0.0;
            A(eq1,iieyp) = 1i*ky/(2.0*h);
            A(eq1,iiezp) = 1i*k/(2.0*h);

            A(eq2,iiexm) = -1i*ky/(2.0*h);
            A(eq2,iieym) = -1.0/(h^2);
            A(eq2,iiezm) = 0.0;
            A(eq2,iiex) = -k0^2*cpdt(2,1,kk);
            A(eq2,iiey) = (k^2 - k0^2*cpdt(2,2,kk)) + 2.0/(h^2);
            A(eq2,iiez) = -ky*k - k0^2*cpdt(2,3,kk);
            A(eq2,iiexp) = 1i*ky/(2.0*h);
            A(eq2,iieyp) = -1.0/(h^2);
            A(eq2,iiezp) = 0.0;

            A(eq3,iiexm) = -1i*k/(2.0*h);
            A(eq3,iieym) = 0.0;
            A(eq3,iiezm) = -1.0/(h^2);
            A(eq3,iiex) = -k0^2*cpdt(3,1,kk);
            A(eq3,iiey) = -ky*k - k0^2*cpdt(3,2,kk);
            A(eq3,iiez) = (ky^2 - k0^2*cpdt(3,3,kk)) + 2.0/(h^2);
            A(eq3,iiexp) = 1i*k/(2.0*h);
            A(eq3,iieyp) = 0.0;
            A(eq3,iiezp) = -1.0/(h^2);
        end

    end

    %--
    % set up rhs vector (current source term)
    rhs = zeros(3*npts,1);
    rhs(1:3:3*npts,1) = 0.;
    rhs(2:3:3*npts,1) = 1i*om*mu0*source';
    rhs(3:3:3*npts,1) = 0.;

    if sparsefill && ~MMS
        
        S = sparse(row,column,v,3*npts,3*npts,(3*npts)*9 - 6);
        
        rf_e = (S)\rhs;
        rf_e = rf_e';
        
    elseif ~sparsefill && ~MMS
        
        %--
        % metallic wall BC
        A(1,1) = 1.0;
        A(2,2) = 1.0;
        A(3,3) = 1.0;
        A(end-2,end-2) = 1.0;
        A(end-1,end-1) = 1.0;
        A(end,end) = 1.0;
        
        rf_e = (A)\rhs;
        rf_e = rf_e';
    
    
    elseif MMS
            
            Ex = 1.0;
            Ey = 1.0;
            Ez = 0.01;
           
            kx = 40.;
            
            ex_solx = Ex*exp(1i*kx*ax);
            ex_soly = Ey*exp(1i*kx*ax);
            ex_solz = Ez*exp(1i*kx*ax);
            
            ex_sol = zeros(1,3*npts);
            ex_sol(1:3:3*npts) = ex_solx;
            ex_sol(2:3:3*npts) = ex_soly;
            ex_sol(3:3:3*npts) = ex_solz;
            
            source_x = zeros(1,npts);
            source_y = zeros(1,npts);
            source_z = zeros(1,npts);
            
        for jj=1:npts
            
            source_x(1,jj) = 1i*ky*(1i*kx*ex_soly(1,jj)) +...
                (ky^2 + kx^2)*ex_solx(1,jj) + 1i*kx*(1i*kx*ex_solz(1,jj)) -...
                k0^2*(cpdt(1,1,jj).*ex_solx(1,jj) + cpdt(1,2,jj).*ex_soly(1,jj) +...
                cpdt(1,3,jj).*ex_solz(1,jj)) + 1i*om*mu0*source(1,jj);
            source_y(1,jj) = 1i*ky*(1i*kx*ex_solx(1,jj)) +...
                (kx^2 - kx^2)*ex_soly(1,jj) - ky*kx*ex_solz(1,jj) -...
                k0^2*(cpdt(2,1,jj).*ex_solx(1,jj) + cpdt(2,2,jj).*ex_soly(1,jj) +...
                cpdt(2,3,jj).*ex_solz(1,jj)) + 1i*om*mu0*source(1,jj);
            source_z(1,jj) = 1i*kx*(1i*kx*ex_solx(1,jj)) +...
                (kx^2 + ky^2)*ex_solz(1,jj) - ky*kx*ex_soly(1,jj) -...
                k0^2*(cpdt(3,1,jj).*ex_solx(1,jj) + cpdt(3,2,jj).*ex_soly(1,jj) +...
                cpdt(3,3,jj).*ex_solz(1,jj)) + 1i*om*mu0*source(1,jj);
            
        end
        
        source_mms = zeros(1,3*npts);
        source_mms(1:3:3*npts) = source_x;
        source_mms(2:3:3*npts) = source_y;
        source_mms(3:3:3*npts) = source_z;
        
        source_mms(1,1) = Ex*exp(1i*kx*axmin); 
        source_mms(1,2) = Ey*exp(1i*kx*axmin);
        source_mms(1,3) = Ez*exp(1i*kx*axmin);
        source_mms(1,end-2) = Ex*exp(1i*kx*axmax); 
        source_mms(1,end-1) = Ey*exp(1i*kx*axmax);
        source_mms(1,end) = Ez*exp(1i*kx*axmax);
        
        rf_e = A\(source_mms' + rhs);
        
        rf_e = rf_e';
        
    else
        
    end

    rf_ex = rf_e(1,1:3:3*npts);
    rf_ey = rf_e(1,2:3:3*npts);
    rf_ez = rf_e(1,3:3:3*npts);

end

