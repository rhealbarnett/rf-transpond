%--------------------------------------------------------------------------------------------------------------%
% solve coupled transport equations                                                                            %
% continuity and momentum eqns                                                                                 %
% NON-CONSERVATIVE FORMS                                                                                           
% (partial derivatives)                                                                                        %
% dvz/dt + vz*d(vz)/dz - nu*d^2(vz)/dz^2 = -(1/mn)(Te+Ti)dn/dz + PF                                                               
% dn/dt + d(n*vz)/dz = Sn                                                                                      %
% staggered n and vz grids                                                                                     %
% momentum eqn upwind convection, central differenced diffusion (IMEX)                                         %
% continuity eqn first order upwind (explicit, flux selecting)                                                 %
% ghost points on density                                                                                      %
% rlbarnett c3149416 140818     
% COMMENTS UPDATED 090719
%--------------------------------------------------------------------------------------------------------------%
%%

% initialise velocity and density 
vx = vx_init;
n = n_init;

%--
% Set 'default' switches.
% sparsefill can be switched to 0, will be slower, but solves. 
% all other switches need to remain as specified. 
staggered = 1;
collocated = 0;
central = 1;
unstable = 0;
sparsefill = 1;

%--
% Set initial BCs. 
rGhost = interp1([nxax(npts-2), nxax(npts-1)], [n(npts-2), n(npts-1)],...
            nxax(npts),'linear','extrap');
lGhost = interp1([nxax(2), nxax(3)], [n(2), n(3)],...
            nxax(1),'linear','extrap');
lvBC_val = LuBC;
rvBC_val = RuBC;

%%
%--------------------------------------------------------------------------------------------------------------%
% INITALISE PLOTS
%--------------------------------------------------------------------------------------------------------------%

if ~MMS && staggered
    vx_source = pressure_source_stag(n_init,const.e,Te,Ti,const.mp,npts,ndx);
    [Ediff, pf] = pond_source({'total',0},{rf_ex,rf_ey,rf_ez},m_s,q_s,'',om,dz,0,{0,''});
    pf_inter = sum(pf,1);
    pf_inter2 = squeeze(sum(pf_inter,2))';
    pf_source = interp1(zax,pf_inter2,vxax,'linear');
    pf_source(1,1) = 0.0; pf_source(1,end) = 0.0;
elseif ~MMS && collocated
    vx_source = pressure_source_col(n_new,const.e,Te,Ti,const.mp,npts-1,dx);
elseif MMS && staggered
    vx_source = mms_source_mom(om,ux,kux,vxax,dt,0,nu,ex_solu,nxax,knx,nx,ex_soln,npts) +...
                pressure_source_stag(n_new,1,0.5,0.5,1,npts,ndx);
    n_source = mms_source_cont(om,nx,knx,nxax(2:npts-1),dt,0,ex_solu(2:npts-1),...
        kux,ux,vxax(2:npts-1),ex_soln(2:npts-1));
    n_source = [0, n_source, 0];
end

if plots
    figure(1)
    set(gcf,'Position',[563 925 560 420])
    plot(nxax(2:npts-1),n_init(2:npts-1),'DisplayName',['time = 0s'])
    xlabel('Position (m)','Fontsize',16)
    ylabel('Density (m^{-3})','Fontsize',16)
    legend('show','Location','west')
    grid on
    hold on

    figure(2)
    set(gcf,'Position',[7 925 560 420])
    if ~MMS
        plot(vxax,vx_init/cs,'DisplayName',['time = 0s'])
    elseif MMS
        plot(vxax,vx_new,'DisplayName',['time = 0s'])
    end
    xlabel('Position (m)','Fontsize',16)
    ylabel('Mach number','Fontsize',16)
    legend('show','Location','southeast')
    grid on
    hold on

    figure(3)
    set(gcf,'Position',[3 476 560 420])
    plot(vxax(2:npts-2),(vx_source(2:npts-2)),'DisplayName',['time = 0s'])
    if MMS
        hold on
        plot(vxax,(pressure_source_stag(n_new,1,nx/2,nx/2,1,npts,ndx)),'DisplayName',['time = 0s'])
    end
    xlabel('Position (m)','Fontsize',16)
    ylabel('Velocity source (ms^{-1})','Fontsize',16)
    legend('show','Location','northwest')
    grid on
    hold on

    figure(4)
    set(gcf,'Position',[563 476 560 420])
    plot(nxax(2:npts-1),n_source(2:npts-1)*dt,'DisplayName',['time = 0s'])
    xlabel('Position (m)','Fontsize',16)
    ylabel('Density source (m^{-3})','Fontsize',16)
    legend('show','Location','northwest')
    grid on
    hold on

    if ~MMS
        figure(5)
        plot(vxax(2:npts-1),pf_source(2:npts-1)*dt,'DisplayName',['time = 0s'])
        xlabel('Position (m)','Fontsize',16)
        ylabel('Ponderomotive source (ms^{-1})','Fontsize',16)
        legend('show','Location','northwest')
        hold on
    end
end
    
%
%--------------------------------------------------------------------------------------------------------------%
% START TIME STEPPING
%--------------------------------------------------------------------------------------------------------------%

vx_rms = zeros(1,nmax);
n_rms = zeros(1,nmax);

vx_I = sparse(eye(npts-1,npts-1));
nI = sparse(eye(npts,npts));

if ~sparsefill
    nA = sparse(npts,npts);
    vx_pos = sparse(npts-1,npts-1);
    vx_neg = sparse(npts-1,npts-1);
    vx_diff = sparse(npts-1,npts-1); 
end

tic

for ii=1:nmax
    
    %--
    % Set exact density and velocity solutions for time dependent
    % MMS calculations.
    if MMS && TD
        ex_solu = u0 + ux*cos(kux*vxax.^2 + om*dt*ii);
        ex_soln = n0 + nx*sin(knx*nxax.^2 + om*dt*ii);
        if continuity && ~momentum
            vx_new = ex_solu;
        end
    end
    
    %--
    % Switch to couple to the wave solver 
    if couple
        
        %--
        % Interpolate density from variable to uniform grid
        n_new_uni = interp1(nxax,n_new,zax,'linear');

        %--
        % Update cold plasma dielectric tensor based on the density
        % calculated from the previous loop
        [om_c,om_p,cpdt,s_arr,d_arr,p_arr] = dielec_tens(q_s,B0,n_new_uni,m_s,om,eps0,npts,1);
        
        %--
        % Ramp the RF electric field amplitude over the first 1000
        % iterations. 
        if ii<=1000
            source_ramp = 1.0/(1001-ii);
            [A,rf_e,rf_ex,rf_ey,rf_ez] = wave_sol(zax,ky,kx,k0,...
            om,mu0,cpdt,source_ramp*source,0,1,1,0);
        else
            [A,rf_e,rf_ex,rf_ey,rf_ez] = wave_sol(zax,ky,kx,k0,...
            om,mu0,cpdt,source,0,1,1,0);
        end
        
        %--
        % Interpolate the electric field from uniform to variable grid
        Ex = interp1(zax,rf_ex,vxax,'linear');
        Ey = interp1(zax,rf_ey,vxax,'linear');
        Ez = interp1(zax,rf_ez,vxax,'linear');
        
    end
    
    %--
    % Start filling the density coefficient matrix for a staggered grid.
    if staggered && (continuity || ~MMS)
        
        %--
        % Fill coefficient matrix using Matlab's sparsefilling... routine?
        if sparsefill
            
            row = zeros(1,(2*npts)-2);
            column = zeros(1,(2*npts)-2);
            n_sparse = zeros(1,(2*npts)-2);
            row(1,1) = 1;
            column(1,1) = 1;
            row(1,end) = npts;
            column(1,end) = npts;
            n_sparse(1,1) = 1;
            n_sparse(1,end) = 1;
        
            for jj=2:npts-1
            
                row(1,2*jj-2) = jj;
                row(1,2*jj-1) = jj;
                column(1,2*jj-2) = jj;

                %--
                % Use the average velocity over the two closest locations
                % to determine whether to use up or down wind scheme. 
                if ((vx(1,jj-1)+vx(1,jj))/2)>0 
                    column(1,2*jj-1) = jj-1;
                    n_sparse(1,2*jj-2) = - (1.0/ndx(1,jj-1))*vx(1,jj);
                    n_sparse(1,2*jj-1) = (1.0/ndx(1,jj-1))*vx(1,jj-1);
                    if MMS
                        n_source(1,jj) = mms_source_cont(om,nx,knx,nxax(1,jj),dt,ii,...
                        ex_solu(1,jj),kux,ux,vxax(1,jj),ex_soln(1,jj));
                    end
                elseif ((vx(1,jj-1)+vx(1,jj))/2)<0
                    column(1,2*jj-1) = jj+1;
                    n_sparse(1,2*jj-2) = (1.0/ndx(1,jj))*vx(1,jj-1);
                    n_sparse(1,2*jj-1) = -(1.0/ndx(1,jj))*vx(1,jj);
                    if MMS
                        n_source(1,jj) = mms_source_cont(om,nx,knx,nxax(1,jj),dt,ii,...
                        ex_solu(1,jj-1),kux,ux,vxax(1,jj-1),ex_soln(1,jj));
                    end     
                end
            end
            
            %--
            % Remove any zeros that weren't filled as this will throw off
            % where the matrix is filled.
            ind_rem = find(column==0);
            column(ind_rem) = [];
            row(ind_rem) = [];
            n_sparse(ind_rem) = [];

            S_nA= sparse(row,column,n_sparse,npts,npts,(2*(npts))-2);

            % Calculate coefficient matrix.
            An_exp = nI + dt*S_nA;
            Anx = -S_nA;
            
        %-- 
        % Fill coefficient matrix using regular loop (slow)
        elseif ~sparsefill
            
            for jj=2:npts-1
            
                if ((vx(1,jj-1)+vx(1,jj))/2)>0 
                    nA(jj,jj) = - (1.0/ndx(1,jj-1))*vx(1,jj);
                    nA(jj,jj-1) = (1.0/ndx(1,jj-1))*vx(1,jj-1);
                    if MMS
                        n_source(1,jj) = mms_source_cont(om,nx,knx,nxax(1,jj),dt,ii,...
                        ex_solu(1,jj-1),kux,ux,vxax(1,jj-1),ex_soln(1,jj));
                    end
                elseif ((vx(1,jj-1)+vx(1,jj))/2)<0
                    nA(jj,jj) = (1.0/ndx(1,jj))*vx(1,jj-1);
                    nA(jj,jj+1) = -(1.0/ndx(1,jj))*vx(1,jj);
                    if MMS
                        n_source(1,jj) = mms_source_cont(om,nx,knx,nxax(1,jj),dt,ii,...
                        ex_solu(1,jj),kux,ux,vxax(1,jj),ex_soln(1,jj));
                    end     
                end
                
            end
            
            An_exp = nI + dt*nA;
            Anx = -nA;
            
            %--
            % Set Dirichlet BCs
            An_exp(1,1) = 1.0;
            An_exp(end,end) = 1.0;
            
        end
        
        %--
        % Set 1's in coefficient matrix for Dirichlet boundary conditions.
        Anx(1,1) = 1.0;
        Anx(end,end) = 1.0;
        
        %--
        % Set ends of density source : to exact solution for steady state
        %                            MMS calculations, or
        %                            : to zero for transport simulations.
        if continuity && SS
            n_source(1,1) = ex_soln(1,1);
            n_source(1,end) = ex_soln(1,end);
        elseif ~MMS || (continuity && TD)
            n_source(1,1) = 0.0; n_source(1,end) = 0.0;
        end
        
        %--
        % Set boundary conditions on n : to exact solution for time
        %                              dependent MMS calculations,  or
        %                              : to the ghost values for transport
        %                              simulations. 
        if continuity && TD
            n(1,1) = ex_soln(1,1);
            n(1,end) = ex_soln(1,end);
        elseif ~MMS
            n(1,1) = lGhost;
            n(1,end) = rGhost;
        end
        
        %--
        % Calculate updated value for density.
        if continuity && SS
            n_new = Anx\n_source';
        elseif (continuity && TD) || ~MMS
            n_new = An_exp*n' + dt*n_source';
        end
        
        %--
        % Transpose solution vector so that it is the correct dimensions
        % for the next loop. 
        n_new = n_new';
        
        %--
        % Save old value of n for rms calculations, set new value of n for
        % velocity calculations. 
        n_tol = n;
        n = n_new;
        
        %--
        % Interpolate the density to the ghost cells
        rGhost = interp1([nxax(npts-2), nxax(npts-1)], [n_new(npts-2), n_new(npts-1)],...
            nxax(npts),'linear','extrap');   
        lGhost = interp1([nxax(2), nxax(3)], [n_new(2), n_new(3)],...
            nxax(1),'linear','extrap');

    % Start filling the density coefficient matrix for a co-located grid.
    % -------------------------------------------------------------------%
    %                                                                    %
    %                   NOT TESTED -- 201022 rlbarnett                   %
    %                   Hasn't been used in a long while.                %
    %                                                                    %
    % -------------------------------------------------------------------%
    elseif collocated
 
        for jj=2:npts-2
            if vx(1,jj)>0
                nA(jj,jj) = -mult*vx(1,jj);
                nA(jj,jj-1) = mult*vx(1,jj-1);
            elseif vx(1,jj)<0
                nA(jj,jj) = mult*vx(1,jj);
                nA(jj,jj+1) = -mult*vx(1,jj+1);
            end
        end
        
        nA(end,end) = -mult*vx(1,end);
        nA(end,end-1) = mult*vx(1,end-1);
        
%         An_exp = nI + dt*nA;
        An_imp = nI - dt*nA;
        
        n(1,1) = lnBC_val;
        
        n_source = n.*n_neut*rate_coeff;

        source_int = trapz(nxax, n_source);

        rflux = vx(end)*n(end);
        ns_mult = rflux/source_int;
        n_source = n_source*ns_mult;
        
        n_source(1,1) = 0.0;
        
        % implicit calculation
        n_new_imp = An_imp\(n' + dt*n_source');
        
        % transpose solution vector
        n_new = n_new_imp;
        n_new = n_new';
    %---------------------------------------------------------------------%
    %                                                                     %
    %                End of untested/update required loop                 %
    %                                                                     %
    %---------------------------------------------------------------------%
    end
    
    %--
    % Fill coefficient matrices for velocity update equation 
    % Positive for v>0 and negative for v<0 on the convective term; 
    % differencing of the diffusion term is central and not 
    % dependent on flow direction
    if momentum || ~MMS
        
        %--
        % Fill coefficient matrix using Matlab's sparsefilling... routine?
        if sparsefill
            
            count = 0;

            row_diff = zeros(1,(3*(npts-1))-4);
            column_diff = zeros(1,(3*(npts-1))-4);
            row_adv = zeros(1,(2*(npts-1))-2);
            column_adv = zeros(1,(2*(npts-1))-2);
            row_diff(1,1) = 1;
            column_diff(1,1) = 1;
            row_diff(1,end) = npts-1;
            column_diff(1,end) = npts-1;
            row_adv(1,1) = 1;
            column_adv(1,1) = 1;
            row_adv(1,end) = npts-1;
            column_adv(1,end) = npts-1;
            vx_sparse_diff = zeros(1,(3*(npts-1))-4);
            vx_sparse_adv = zeros(1,(2*(npts-1))-2);
            
            for jj=2:npts-2            
            
                row_diff(1,jj+count*2) = jj;
                row_diff(1,jj+count*2+1) = jj;
                row_diff(1,jj+count*2+2) = jj;
                column_diff(1,jj+2*count) = jj;
                column_diff(1,jj+2*count+1) = jj-1;
                column_diff(1,jj+2*count+2) = jj+1;

                row_adv(1,2*jj-2) = jj;
                row_adv(1,2*jj-1) = jj;
                column_adv(1,2*jj-2) = jj;

                %--
                % Fill the diffusion coefficient matrix.
                vx_sparse_diff(1,jj+2*count) = - (1.0/(vdx(1,jj-1)*vdx(1,jj)))*(2.0*nu);
                vx_sparse_diff(1,jj+2*count+1) = (2.0/(vdx(1,jj-1)*(vdx(1,jj) + vdx(1,jj-1))))*nu;
                vx_sparse_diff(1,jj+2*count+2) = (2.0/((vdx(1,jj-1) + vdx(1,jj))*vdx(1,jj)))*nu;

                count = count+1;

                %--
                % Check for sign of velocity to fill the advection
                % coefficient terms.
                if vx(1,jj)>0
                    column_adv(1,2*jj-1) = jj-1;
                    if ~MMS
                        vx_sparse_adv(1,2*jj-2) = - (1.0/vdx(1,jj-1))*vx(1,jj) -...
                            (1.0/n(1,jj))*n_source(1,jj);
                    elseif MMS
                        vx_sparse_adv(1,2*jj-2) = - (1.0/vdx(1,jj-1))*vx(1,jj);
                    end
                    vx_sparse_adv(1,2*jj-1) = (1.0/vdx(1,jj-1))*vx(1,jj);
                elseif vx(1,jj)<0
                    column_adv(1,2*jj-1) = jj+1;
                    if ~MMS
                        vx_sparse_adv(1,2*jj-2) = (1.0/vdx(1,jj))*vx(1,jj) -...
                            (1.0/n(1,jj+1))*n_source(1,jj+1);
                    elseif MMS
                        vx_sparse_adv(1,2*jj-2) = (1.0/vdx(1,jj))*vx(1,jj);
                    end
                    vx_sparse_adv(1,2*jj-1) = - (1.0/vdx(1,jj))*vx(1,jj);
                end
                    
            end
            
            %--
            % Remove any zero indices as this will skew the coefficient
            % matrix.
            ind_rem = find(column_adv==0);
            column_adv(ind_rem) = [];
            row_adv(ind_rem) = [];
            vx_sparse_adv(ind_rem) = [];

            S_vxE = sparse(row_adv,column_adv,vx_sparse_adv,npts-1,npts-1,(2*(npts-1))-2);
            S_vxI = sparse(row_diff,column_diff,vx_sparse_diff,npts-1,npts-1,(3*(npts-1))-4);  

            %--
            % Calculate the advection (Avx_exp) and diffusion (Avx_imp)
            % matrices.
            Avx_exp = vx_I + dt*S_vxE;
            Avx_imp = vx_I - dt*S_vxI;
            
            %--
            % Avx is used for steady state MMS testing.
            Avx = -(S_vxE + S_vxI);
         
        %-- 
        % Fill coefficient matrix using regular loop (slow).    
        elseif ~sparsefill
            
            for jj=2:npts-2 
                
                if vx(1,jj)>0
                    vx_pos(jj,jj) = - (1.0/vdx(1,jj-1))*vx(1,jj) -...
                        (1.0/n(1,jj))*n_source(1,jj);
                    vx_pos(jj,jj-1) = (1.0/vdx(1,jj-1))*vx(1,jj);
                elseif vx(1,jj)<0
                    vx_neg(jj,jj) = (1.0/vdx(1,jj))*vx(1,jj) -...
                        (1.0/n(1,jj+1))*n_source(1,jj+1);
                    vx_neg(jj,jj+1) = - (1.0/vdx(1,jj))*vx(1,jj);
                end
                vx_diff(jj,jj) = - (2.0/(vdx(1,jj-1)*vdx(1,jj)))*(nu);
                vx_diff(jj,jj-1) = (2.0/(vdx(1,jj-1)*(vdx(1,jj) + vdx(1,jj-1))))*nu;
                vx_diff(jj,jj+1) = (2.0/((vdx(1,jj-1) + vdx(1,jj))*vdx(1,jj)))*nu;
                
            end
                
            %--
            % Can switch between advection only and advection-diffusion
            % equation. 
            if upwind 
                vxA = vx_pos + vx_neg;
            elseif central
                vxE = vx_pos + vx_neg;
                vxI = vx_diff;
            end
            
            %--
            % Calculate advection (Avx_exp) and diffusion (Avx_imp)
            % coefficient matrices. 
            Avx_exp = vx_I + dt*vxE;
            Avx_imp = vx_I - dt*vxI;
            
            %--
            % Avx is used for steady state MMS testing.
            Avx = -(vxE + vxI);
            
            %--
            % Set coefficients for Dirichlet boundary conditions.
            Avx_exp(1,1) = 1.0; Avx_exp(end,end) = 1.0;
            Avx_imp(1,1) = 1.0; Avx_imp(end,end) = 1.0;

        end
        
        %--
        % Set Dirichlet boundary conditions for steady state MMS. 
        Avx(1,1) = 1.0;
        Avx(end,end) = 1.0;

        %--
        % Set Dirichlet boundary conditions based on either MMS (exact
        % solution at boundaries) or transport simulations (l/rvBC_val).
        if momentum && TD
            vx(1,1) = ex_solu(1,1);
            vx(1,end) = ex_solu(1,end);
        elseif ~MMS
            vx(1,1) = lvBC_val;
            vx(1,end) = rvBC_val;
        end

        %--
        % Calculate the right hand side 'source' terms for a staggered
        % grid.
        if staggered && ~MMS
            %--
            % Grad-P source, and set boundaries to zero (BCs set in
            % coefficient matrix and vx)
            vx_source = pressure_source_stag(n,const.e,Te,Ti,const.mp,npts,ndx);
            vx_source(1,1) = 0.0; vx_source(1,end) = 0.0;
            
            %--
            % Calculate the ponderomotive source term.
            if test && ii<=1000
                source_ramp = 1.0/(1001-ii);
                [Ediff, pf] = pond_source({'para',0},{rf_ex*source_ramp,...
                    rf_ey*source_ramp,rf_ez*source_ramp},m_s,q_s,'',om,dz,0,{0,''});
            else
                [Ediff, pf] = pond_source({'total',0},{rf_ex,rf_ey,rf_ez},m_s,q_s,'',om,dz,0,{0,''});
            end
            pf_inter = sum(pf,1);
            pf_inter2 = squeeze(sum(pf_inter,2))';
            pf_source = interp1(zax,pf_inter2,vxax,'linear');
            
            %--
            % Set boundaries of source to zero (velocity BCs set in
            % coefficient matrix and vx)
            pf_source(1,1) = 0.0; pf_source(1,end) = 0.0;
        
        %--
        % Calculate MMS source term. 
        elseif staggered && momentum
            vx_source = mms_source_mom(om,ux,kux,vxax,dt,ii,nu,ex_solu,nxax,knx,nx,ex_soln,npts) +...
                pressure_source_stag(n_new,1,0.5,0.5,1,npts,ndx);
            %--
            % Set MMS boundary conditions depending on whether it is a
            % steady state (SS) or time dependent (TD) test.
            if SS
                vx_source(1,1) = ex_solu(1,1);
                vx_source(1,end) = ex_solu(1,end);
            elseif TD
                vx_source(1,1) = 0;
                vx_source(1,end) = 0;
            end
            
        %--
        % Calculate source term for a co-located grid 
        % CURRENTLY OBSOLETE 201027 rlbarnett :)
        elseif collocated
            vx_source = pressure_source_col(n,const.e,Te,Ti,m,npts-1,dx);
        end

        %--
        % Calculate solution for MMS steady state, MMS time dependent or
        % transport simulation. 
        if momentum && SS
            vx_new = Avx\vx_source';
        elseif (momentum && TD)
            vx_newE = Avx_exp*vx';
            vx_new = Avx_imp\(vx_newE + dt*(vx_source'));
        elseif ~MMS 
            vx_newE = Avx_exp*vx';
            vx_new = Avx_imp\(vx_newE + dt*(vx_source' - pf_source'));
        end
        
        %--
        % Transpose solution vector so that it is the correct shape for the
        % next loop.
        vx_new = vx_new';
        
        
        %--
        % MMS steady state check for convergence. 
        if MMS && SS
            if rms(vx_new - vx)<=tol && rms(n_new - n_tol)<=tol
                fprintf('tolerance reached, ii=%d\n',ii)
                fprintf('velocity rms error = %d\n', rms(vx - vx_new))
                fprintf('density rms error = %d\n', rms(n_tol - n_new))
                return
            elseif mod(ii,round(nmax/10))==0 || ii==nmax
                fprintf('velocity rms error = %d\n', rms(vx - vx_new))
                fprintf('density rms error = %d\n', rms(n_tol - n_new))
            end
        end
        
        vx = vx_new;
        
        if test && rms(vx_new - vx)/rms(vx)<=1.0e-8 && rms(n_new - n_tol)/rms(n_tol)<=1.0e-8
            fprintf('tolerance reached, ii=%d\n',ii)
            return
        end
        
    end
    
    %--
    % Check for NaNs in the velocity array (can be density, whichever) and
    % stop iterations if any exist.
    nan_check = isnan(vx_new);
    
    if sum(nan_check) ~= 0
        unstable = 1;
        fprintf('unstable, ii=%d\n',ii)
        return
    end

    %--
    % Plot all fields if plot switch is on. 
    if plots && (mod(ii,round(nmax/5))==0 || ii==nmax)
        
        %--
        % Plot density (include exact solution for MMS)
        figure(1)
        set(gcf,'Position',[563 925 560 420])
        plot(nxax,n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax+ndx(1)) max(nxax-ndx(end))])
        hold on
        if MMS
            plot(nxax,ex_soln,'--','DisplayName',['exact = ' num2str(double(ii)*dt) ' s'])
        end
        
        %--
        % Plot mach number (include exact solution for MMS)
        figure(2)
        set(gcf,'Position',[7 925 560 420])
        plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        if MMS
            plot(vxax,vx_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
            plot(vxax,ex_solu,'--','DisplayName',['exact = ' num2str(double(ii)*dt) ' s'])
        end
        xlim([min(vxax) max(vxax)])
        hold on
        
        %--
        % Plot the velocity source (grad-p only).
        figure(3)
        set(gcf,'Position',[3 476 560 420])
        plot(vxax(2:npts-2),(vx_source(2:npts-2)),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        if MMS
            hold on
            plot(vxax,(pressure_source_stag(n_new,1,nx/2,nx/2,1,npts,ndx)),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        end
        xlim([min(vxax) max(vxax)])
        hold on
        
        %--
        % Plot the density source (should remain constant over each loop).
        figure(4)
        set(gcf,'Position',[563 476 560 420])
        plot(nxax(2:npts-1),n_source(2:npts-1)*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlabel('Position (m)','Fontsize',16)
        ylabel('Density source ms^{-1}','Fontsize',16)
        legend('show','Location','northwest')
        hold on
        
        
        if ~MMS
            %-- 
            % Plot ponderomotive source term.
            figure(5)
            plot(vxax(2:npts-1),pf_source(2:npts-1)*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
            xlabel('Position (m)','Fontsize',16)
            ylabel('Ponderomotive source (ms^{-1})','Fontsize',16)
            legend('show','Location','northwest')
            hold on
        end
        
        
    elseif (mod(ii,save_iter)==0 && sfile) || unstable || (mod(ii,nmax)==0 && sfile)
        
        transport.dt = dt;
        transport.n_source = n_source;
        transport.vx_source = vx_source;
        transport.pf_source = pf_source;
        transport.vx_new = vx_new;
        transport.n_new = n_new; 
        transport.ii = ii;
        transport.cs = cs;
        transport.vxax = vxax;
        transport.nxax = nxax;
        transport.vdx = vdx;
        transport.ndx = ndx;
        transport.Te = Te;
        transport.Ti = Ti;
        transport.npts = npts;
        transport.tmax = tmax;
        transport.xmin = zmin;
        transport.xmax = zmax;
        transport.nu = nu;
        transport.freq = freq;
        transport.E_0 = E_0;
        transport.rf_ex = rf_ex;
	    transport.rf_ey = rf_ey;
	    transport.rf_ez = rf_ez;
	    transport.zax = zax;
        transport.pond = pf;
        transport.pond_summed = pf_source;
        if ~test
            transport.poyn = poyn;
            transport.source = source;
            transport.kx = kx;
            transport.ky = ky;
            transport.period = period;
            transport.B0 = B0;
        end

        filename = strcat(filepath,'outputs/coupled_',num2str(ii),'.mat');
        save(filename,'-struct','transport');
        
        continue
    end    
end

toc

