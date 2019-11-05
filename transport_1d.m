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

% import parameter file
% params_transport_wave_ACM;
% transport_vardx;
% transport_test;
% transport_mms;
lapd_equib;

% initialise velocity and density 
vx = vx_new;
n = n_new;

% staggered = NaN;
% collocated = NaN;
% v_ldirichlet = NaN;
% v_rdirichlet = NaN;
% v_rneumann = NaN;
% v_lneumann = NaN;
% v_periodic = NaN;
% n_ldirichlet = NaN;
% n_rdirichlet = NaN;
% n_rneumann = NaN;
% n_lneumann = NaN;
% n_periodic = NaN;
% explicit = NaN;
% implicit = NaN;
% MMS = NaN;
% momentum = NaN;
% continuity = NaN;

staggered = 1;
collocated = 0;
v_ldirichlet = 1;
v_rdirichlet = 1;
v_rneumann = 0;
v_lneumann = 0;
v_periodic = 0;
n_ldirichlet = 0;
n_rdirichlet = 0;
n_rneumann = 0;
n_lneumann = 0;
n_periodic = 0;
MMS = 0;
momentum = 0;
continuity = 0;
central = 1;
upwind = 0;
unstable = 0;
plots = 1;


rGhost = interp1([nxax(npts-2), nxax(npts-1)], [n_new(npts-2), n_new(npts-1)],...
            nxax(npts),'linear','extrap');
lGhost = interp1([nxax(2), nxax(3)], [n_new(2), n_new(3)],...
            nxax(1),'linear','extrap');
lvBC_val = LuBC;
rvBC_val = RuBC;

% test_type = 'run MMS (yes/no)? ';
% testt = input(test_type, 's');
% if isempty(testt)
%     testt = 'no';
% end
% 
% if strcmp(testt,'yes')
%     MMS = 1;
%     mms_type = 'continuity, momentum, or coupled? ';
%     mmst = input(mms_type, 's');
%     if strcmp(mmst,'continuity')
%         continuity = 1;
%         momentum = 0;
%     elseif strcmp(mmst,'momentum')
%         continuity = 0;
%         momentum = 1;
%     elseif strcmp(mmst,'coupled')
%         continuity = 1;
%         momentum = 1;
%     end
% elseif strcmp(testt,'no')
%     MMS = 0;
%     continuity = 0;
%     momentum = 0;
%     coupled = 0;
% end
% 
% grid_type = 'staggered or collocated grid? ';
% gridt = input(grid_type, 's');
% if isempty(gridt)
%     gridt = 'staggered';
% end
% 
% if strcmp(gridt,'staggered')
%     staggered = 1;
%     collocated = 0;
% elseif strcmp(gridt,'collocated')
%     staggered = 0;
%     collocated = 1;
% end
% 
% if (isnan(staggered)) || (isnan(collocated))
%     error("Check spelling and/or type of answer for %s.\n",'"staggered or collocated grid?"')
%     return
% end
% 
% if staggered
%     
%     ln_bound_type = 'Left (ghost) BC type? (dirichlet, neumann, periodic, linear extrap) ';
%     leftGhost = input(ln_bound_type, 's');
%     if isempty(leftGhost)
%         leftGhost = 'linear extrap';
%     end
%     
%     if strcmp('periodic',leftGhost)
% 
%     else
%         rn_bound_type = 'Right (ghost) BC type? (dirichlet, neumann, linear extrap) ';
%         rightGhost = input(rn_bound_type, 's');
%         if isempty(rightGhost)
%             rightGhost = 'linear extrap';
%         end
%     end
%     
%     if strcmp('linear extrap',leftGhost)
%         lGhost = interp1([nxax(2), nxax(3)], [n_new(2), n_new(3)],...
%             nxax(1),'linear','extrap');
%         n_ldirichlet = 0;              
%         n_lneumann = 0;
%         n_periodic = 0;
%     elseif strcmp('dirichlet',leftGhost)
%         n_ldirichlet = 1;              
%         n_lneumann = 0;
%         n_periodic = 0;
%     elseif strcmp('neumann',leftGhost)
%         n_ldirichlet = 0;                
%         n_lneumann = 1; 
%         n_periodic = 0;
%     end  
%     
%     if (isnan(n_ldirichlet)) || (isnan(n_lneumann)) || (isnan(n_periodic))
%         error("Check spelling and/or type of answer for %s.\n",...
%             '"Left (ghost) BC type? (dirichlet, neumann, periodic, linear extrap)"')
%         return
%     end
%     
%     if strcmp('linear extrap',leftGhost)
%     else
%         ln_bound_val = 'Left (ghost) BC value for density? ';
%         lnBC_val = input(ln_bound_val);
%         if isempty(lnBC_val)
%             lnBC_val = LnBC;
%         end
%     end
%     
%     if strcmp('linear extrap',rightGhost)
%         rGhost = interp1([nxax(npts-2), nxax(npts-1)], [n_new(npts-2), n_new(npts-1)],...
%             nxax(npts),'linear','extrap');
%         n_rdirichlet = 0;                
%         n_rneumann = 0; 
%         n_periodic = 0;
%     elseif strcmp('dirichlet',rightGhost)
%         n_rdirichlet = 1;              
%         n_rneumann = 0;
%         n_periodic = 0;
%     elseif strcmp('neumann',rightGhost)
%         n_rdirichlet = 0;                
%         n_rneumann = 1; 
%         n_periodic = 0;
%     end
%     
%     if strcmp('linear extrap',rightGhost)
%     else
%         rn_bound_val = 'Right (ghost) BC value for density? ';
%         rnBC_val = input(rn_bound_val);
%         if isempty(rnBC_val)
%             rnBC_val = RnBC;
%         end
%     end
%     
%     if (isnan(n_rdirichlet)) || (isnan(n_rneumann)) || (isnan(n_periodic))
%         error("Check spelling and/or type of answer for %s.\n",...
%             '"Right (ghost) BC type? (dirichlet, neumann, periodic)"')
%         return
%     end
%     
% end
% 
% if collocated
%     
%     if vx_new(1,2) > 0
%         promptlnBC = 'Left BC type for density? (dirichlet, neumann, periodic) ';
%         leftnBC = input(promptlnBC, 's');
%         if strcmp('dirichlet',leftnBC)
%             n_ldirichlet = 1;              
%             n_lneumann = 0;
%             n_periodic = 0;
%         elseif strcmp('neumann',leftnBC)
%             n_ldirichlet = 0;                
%             n_lneumann = 1; 
%             n_periodic = 0;
%         end
%         promptlnBCval = 'Left BC value for density? ';
%         lnBC_val = input(promptlnBCval);
%     elseif vx_new(1,2) < 0 
%         fprintf("Left BC not required on density for the given flux direction.\n")
%         n_ldirichlet = 0;
%         n_lneumann = 0;
%         n_periodic = 0;
%     end
%     if vx_new(1,end-1) < 0
%         promptrnBC = 'Right BC type for density? (dirichlet, neumann, periodic) ';
%         rightnBC = input(promptrnBC, 's');
%         if strcmp('dirichlet',rightnBC)
%             n_rdirichlet = 1;              
%             n_rneumann = 0;
%             n_periodic = 0;
%         elseif strcmp('neumann',rightnBC)
%             n_rdirichlet = 0;                
%             n_rneumann = 1; 
%             n_periodic = 0;
%         end
%         rn_bound_val = 'Right BC value for density? ';
%         rnBC_val = input(rn_bound_val);
%     elseif vx_new(1,end-1) > 0 
%         fprintf("Right BC not required on density for the given flux direction.\n")
%         n_rdirichlet = 0;
%         n_rneumann = 0;
%         n_periodic = 0;
%     end
% end
% 
% %%
% 
% upwind = NaN;
% central = NaN;
% 
% spatial_scheme = 'Spatial differencing scheme for momentum equation? (upwind or central) ';
% scheme = input(spatial_scheme, 's');
% if isempty(scheme)
%     scheme = 'central';
% end
% 
% if strcmp(scheme,'upwind')
%     upwind = 1;
%     central = 0;
% elseif strcmp(scheme, 'central')
%     upwind = 0;
%     central = 1;
% end
% 
% % if central && nu==0
% %     fprintf("nu==0: central difference scheme not stable.\n")
% %     return
% % end
% 
% %%
% %--------------------------------------------------------------------------------------------------------------%
% % select boundary conditions 
% 
% 
% %----- need a prompt here to check whether the velocity is positive or
% % negative next to the boundary, as this will determine whether left or right (or both)
% % BC is required
% 
% if central
% 
%     lv_bound_type = 'Left BC type for velocity? (dirichlet, neumann, periodic) ';
%     leftvBC = input(lv_bound_type, 's');
%     if isempty(leftvBC)
%         leftvBC = 'dirichlet';
%     end
%     if strcmp('periodic',leftvBC)
% 
%     else
%         rv_bound_type = 'Right BC type for velocity? (dirichlet or neumann) ';
%         rightvBC = input(rv_bound_type, 's');
%         if isempty(rightvBC)
%             rightvBC = 'dirichlet';
%         end
%     end
% 
%     if strcmp('dirichlet',leftvBC)
%         v_ldirichlet = 1;               % this allows the use of ldirichlet as a logical 
%         v_lneumann = 0;
%         v_periodic = 0;
%     elseif strcmp('neumann',leftvBC)
%         v_ldirichlet = 0;                
%         v_lneumann = 1; 
%         v_periodic = 0;
%     end
% 
%     if strcmp('dirichlet',rightvBC)
%         v_rdirichlet = 1;
%         v_rneumann = 0;
%         v_periodic = 0;
%     elseif strcmp('neumann',rightvBC)
%         v_rneumann = 1;
%         v_rdirichlet = 0;
%         v_periodic = 0;
%     end
% 
%     if strcmp('periodic',leftvBC)
%         v_rdirichlet = 0;
%         v_rneumann = 0;
%         v_ldirichlet = 0;
%         v_lneumann = 0;
%         v_periodic = 1;
%     end
% 
%     lv_bound_val = 'Left BC value for velocity? ';
%     lvBC_val = input(lv_bound_val);
%     if isempty(lvBC_val)
%         lvBC_val = LuBC;
%     end
%     if strcmp('periodic',leftvBC)
% 
%     else
%         rv_bound_val = 'Right BC value for velocity? ';
%         rvBC_val = input(rv_bound_val);
%         if isempty(rvBC_val)
%             rvBC_val = RuBC;
%         end
%     end
% 
%     if v_ldirichlet && v_rdirichlet
%         fprintf("left and right velocity BCs are dirichlet; lBC = %f, rBC = %f.\n", lvBC_val, rvBC_val)
%     elseif v_lneumann && v_rdirichlet
%         fprintf("left velocity BC is neumann, lBC = %f; right velocity BC is dirichlet, rBC = %f.\n", lvBC_val, rvBC_val)
%     elseif v_rneumann && v_ldirichlet
%         fprintf("left velocity BC is dirichlet, lBC = %f; right velocity BC is neumann, rBC = %f.\n", lvBC_val, rvBC_val)
%     elseif v_rneumann && v_lneumann
%         fprintf("left and right velocity BCs are dirichlet; lBC = %f, rBC = %f.\n", lvBC_val, rvBC_val)
%     elseif v_periodic
%         fprintf("periodic velocity BCs\n")
%     end
%     
% elseif upwind
%     
%     if vx_new(1,2) > 0
%             lv_bound_type = 'Left BC type for velocity? (dirichlet, neumann, periodic) ';
%             leftvBC = input(lv_bound_type, 's');
%             if strcmp('dirichlet',leftvBC)
%                 v_ldirichlet = 1;              
%                 v_lneumann = 0;
%                 v_periodic = 0;
%             elseif strcmp('neumann',leftvBC)
%                 v_ldirichlet = 0;                
%                 v_lneumann = 1; 
%                 v_periodic = 0;
%             end
%             lv_bound_val = 'Left BC value for velocity? ';
%             lvBC_val = input(lv_bound_val);
%     elseif vx_new(1,2) < 0 
%             fprintf("Left BC not required\n")
%             v_ldirichlet = 0;
%             v_lneumann = 0;
%             v_periodic = 0;
%     end
%     
%     if vx_new(1,end-1) < 0
%             rv_bound_type = 'Right BC type? (dirichlet, neumann, periodic) ';
%             rightvBC = input(rv_bound_type, 's');
%             if strcmp('dirichlet',rightvBC)
%                 v_rdirichlet = 1;              
%                 v_rneumann = 0;
%                 v_periodic = 0;
%             elseif strcmp('neumann',rightvBC)
%                 v_rdirichlet = 0;                
%                 v_rneumann = 1; 
%                 v_periodic = 0;
%             end
%             rv_bound_val = 'Right BC value? ';
%             rvBC_val = input(rv_bound_val);
%     elseif vx_new(1,end-1) > 0 
%             fprintf("Right BC not required\n")
%             v_rdirichlet = 0;
%             v_rneumann = 0;
%             v_periodic = 0;
%     end
%     
% end

%%
%--------------------------------------------------------------------------------------------------------------%
% INITALISE PLOTS; INCLUDE INITIAL CONDITIONS
%--------------------------------------------------------------------------------------------------------------%

if ~MMS && staggered
    vx_source = source_stag(n_new,const.e,Te,Ti,const.mp,npts,ndx);
    pf_source = pond_source(const.me,m,om,const.e,Efield,vdx);
    pf_source = [0,pf_source,0];
elseif ~MMS && collocated
    vx_source = source_col(n_new,const.e,Te,Ti,const.mp,npts-1,dx);
elseif MMS && staggered
    vx_source = mms_source_mom(om,ux,kux,vxax,dt,0,nu,vx_new,nxax,knx,nx,n_new,npts) +...
            source_stag(n_new,1,nx/2,nx/2,1,npts,ndx);
    n_source = mms_source_cont(om,nx,knx,nxax(2:npts-1),dt,0,vx_new(2:npts-1),...
        kux,ux,vxax(2:npts-1),n_new(2:npts-1));
    n_source = [0, n_source, 0];
end

figure(1)
set(gcf,'Position',[563 925 560 420])
% semilogy(nxax(2:npts-1),n_new(2:npts-1),'DisplayName',['time = 0s'])
plot(nxax(2:npts-1),n_new(2:npts-1),'DisplayName',['time = 0s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density (m^{-3})','Fontsize',16)
legend('show','Location','west')
grid on
hold on

figure(2)
set(gcf,'Position',[7 925 560 420])
if ~MMS
    plot(vxax,vx_new/cs,'DisplayName',['time = 0s'])
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
    plot(vxax,(source_stag(n_new,1,nx/2,nx/2,1,npts,ndx)),'DisplayName',['time = 0s'])
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

    figure(6)
    subplot(2,1,1)
    plot(xax,real(rf_ez),'k','DisplayName',['Re[E_{||}], time = 0s'])
    set(gca, 'XTickLabel', [])
    ylabel('E_{||} (Vm^{-1})','Fontsize',16)
    hold on
    plot(xax,imag(rf_ez),'--r','DisplayName',['Im[E_{||}], time = 0s'])
    set(gca, 'XTickLabel', [])
    legend('show','Location','northwest')
    set(gca,'Fontsize',30)
    hold on
    subplot(2,1,2)
    plot(xax,abs(rf_ez).^2,'DisplayName',['time = 0s'])
    xlabel('Position (m)','Fontsize',16)
    ylabel('|E_{||}|^2 (V^2m^{-2})','Fontsize',16)
    legend('show','Location','northwest')
    set(gca,'Fontsize',30)
    hold on
end

%
%--------------------------------------------------------------------------------------------------------------%
% START TIME STEPPING
%--------------------------------------------------------------------------------------------------------------%

count = 2;
timerVal = tic;

vx_rms = zeros(1,nmax);
n_rms = zeros(1,nmax);

for ii=1:nmax
    
    if MMS
        ex_solu = u0 + ux*cos(kux*vxax.^2 + om*dt*ii);
        ex_soln = n0 + nx*sin(knx*nxax.^2 + om*dt*ii);
    end
    
%     set the vectors with the old value going into the next loop
    vx = vx_new;
    rGhost = interp1([nxax(npts-2), nxax(npts-1)], [n_new(npts-2), n_new(npts-1)],...
        nxax(npts),'linear','extrap');   
    lGhost = interp1([nxax(2), nxax(3)], [n_new(2), n_new(3)],...
        nxax(1),'linear','extrap');
    
    n_new_uni = interp1(nxax,n_new,xax,'linear');

    [om_c,om_p,cpdt,s_arr,d_arr,p_arr] = dielec_tens(q_s,B0,n_new_uni,m_s,om,eps0,npts);
    [A,source,rf_e,rf_ex,rf_ey,rf_ez,diss_pow] = wave_sol(xax,ky,kx,k0,...
    om,mu0,cpdt,source_width,source_loc,0,source_mult,1);

    Efield = abs(rf_ez);
    Efield = Efield.^2;
    
    Efield = interp1(xax,Efield,vxax,'linear');
    
    if staggered && (continuity || ~MMS)
        
        % fill n coefficient matrix using the averaged value of the
        % velocities on adjacent grid points
        for jj=2:npts-1
            if ((vx(1,jj-1)+vx(1,jj))/2)>0 
                nA(jj,jj) = - (1.0/ndx(1,jj-1))*vx(1,jj);
                nA(jj,jj-1) = (1.0/ndx(1,jj-1))*vx(1,jj-1);
                if MMS
                    n_source(1,jj) = mms_source_cont(om,nx,knx,nxax(1,jj),dt,ii,...
                    ex_solu(1,jj),kux,ux,vxax(1,jj),ex_soln(1,jj));
                end
            elseif ((vx(1,jj-1)+vx(1,jj))/2)<0
                nA(jj,jj) = (1.0/ndx(1,jj))*vx(1,jj-1);
                nA(jj,jj+1) = -(1.0/ndx(1,jj))*vx(1,jj);
                if MMS
                    n_source(1,jj) = mms_source_cont(om,nx,knx,nxax(1,jj),dt,ii,...
                    ex_solu(1,jj-1),kux,ux,vxax(1,jj-1),ex_soln(1,jj));
                end                    
            end
        end
        
        % set source density ghost points to zero 
        n_source(1,1) = 0.0; n_source(1,end) = 0.0;

        % build full coefficient matrices
        An_exp = nI + dt*nA;
%         An_imp = nI - dt*nA;
        
        % override values in top and bottom rows to reflect neumann
        % boundary conditions for the implicit calculation
        
%         An_imp(1,1) = 1.0; %An_imp(1,2) = -1.0;
%         An_imp(end,end) = 1.0; %An_imp(end,end-1) = -1.0;   
        An_exp(1,1) = 1.0;
        An_exp(end,end) = 1.0;
        
        % calculate explicit solution
%         n_new_exp = An_exp*n' + dt*n_source';
        % directly override solution vector to include neumann boundary
        % conditions for explicit method
%         n_new_exp(1,1) = n_new_exp(2,1);
%         n_new_exp(end,1) = n_new_exp(end-1,1);
        
        % zero old rhs values for top and bottom boundary equations for
        % implicit calculation
        if continuity
            n(1,1) = n0 + nx*sin(knx*min(nxax)^2 + om*dt*ii);
            n(1,end) = n0 + nx*sin(knx*max(nxax)^2 + om*dt*ii);
        elseif ~MMS
            n(1,1) = lGhost;
            n(1,end) = rGhost;
        end
        
        % implicit calculation
%         n_new_imp = An_imp\(n' + dt*n_source');
        n_new = An_exp*n' + dt*n_source';
        
        
        % transpose solution vector
%         n_new = n_new_imp;
        n_new = n_new';
        n = n_new;
        
     
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
        
    end
    
    if momentum || ~MMS
        % fill coefficient matrices for momentum equation, positive for v>0 and
        % negative for v<0 on the convective term; differencing of the
        % diffusion term is central and not dependent on flow direction
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
            vx_diff(jj,jj) = - (1.0/(vdx(1,jj-1)*vdx(1,jj)))*(2.0*nu);
            vx_diff(jj,jj-1) = (2.0/(vdx(1,jj-1)*(vdx(1,jj) + vdx(1,jj-1))))*nu;
            vx_diff(jj,jj+1) = (2.0/((vdx(1,jj-1) + vdx(1,jj))*vdx(1,jj)))*nu;

        end

        % construct full coefficient matrix for momentum equation
        if upwind 
            vxA = vx_pos + vx_neg;
        elseif central
            vxE = vx_pos + vx_neg;
            vxI = vx_diff;
        end

        % build full coefficient matrices
        Avx_exp = vx_I + dt*vxE;
        Avx_imp = vx_I - dt*vxI;
%         Avx_imp = -vxA;
        % override top and bottom rows to include dirichlet boundary conditions
        % for the momentum equation (explicit and implicit methods)
        Avx_exp(1,1) = 1.0; Avx_exp(end,end) = 1.0;
        Avx_imp(1,1) = 1.0; Avx_imp(end,end) = 1.0;

        % override values in top and bottom rows to reflect neumann
        % boundary conditions for the implicit calculation
    %     Avx_imp(1,2) = -1.0;% Avx_imp(end,end-1) = -1.0;
        % ensure that the velocity value at the boundaries is correct
        if momentum
            vx(1,1) = u0 + ux*cos(kux*min(vxax)^2 + om*dt*ii);
            vx(1,end) = u0 + ux*cos(kux*max(vxax)^2 + om*dt*ii);
        elseif ~MMS
            vx(1,1) = lvBC_val;
            vx(1,end) = rvBC_val;
        end

        % calculate the source term
        if staggered && ~MMS
            vx_source = source_stag(n,const.e,Te,Ti,const.mp,npts,ndx);
            pf_source = pond_source(const.me,m,om,const.e,Efield,vdx);
            pf_source = [0,pf_source,0];
        elseif staggered && momentum
            vx_source = mms_source_mom(om,ux,kux,vxax,dt,ii,nu,ex_solu,nxax,knx,nx,ex_soln,npts) +...
                source_stag(n,1,nx/2,nx/2,1,npts,ndx);
            vx_source(1,1) = u0 + ux*cos(kux*min(vxax)^2 + om*dt*ii);
            vx_source(1,end) = u0 + ux*cos(kux*max(vxax)^2 + om*dt*ii);
        elseif collocated
            vx_source = source_col(n,const.e,Te,Ti,m,npts-1,dx);
        end

        % zero the source term at the boundaries as it is not used (dirichlet
        % boundary conditions will override the source)
        vx_source(1,1) = 0.0;
        vx_source(1,end) = 0.0;

        % explicit calculation
    %     vx_new_exp = Avx_exp*vx' + dt*(vx_source' + pf_source');
        % implicit calculation
%         vx_new_imp = Avx_imp\(vx' + dt*(vx_source' - pf_source'));
%         vx_new_imp = Avx_imp\(vx' + dt*vx_source');
%         vx_new_imp = Avx_imp\(vx_source');
        vx_newE = Avx_exp*vx';
        vx_new = Avx_imp\(vx_newE + dt*(vx_source' + pf_source'));

        % transpose solution vector
%         vx_new = vx_new_imp;
        vx_new = vx_new';
    end
    
%     reset CFL condition based on the lowest dt out of the
%     convective/diffusive CFLs
%     if (cfl_fact*(dx^2)/(2.0*nu))<(cfl_fact*dx/max(abs(vx_new)))
%         dt = cfl_fact*(dx^2)/(2.0*nu);
%     elseif (cfl_fact*(dx^2)/(2.0*nu))>(cfl_fact*dx/max(abs(vx_new)))
%         dt = cfl_fact*dx/max(abs(vx_new));
%     end

%     will stop running script if either of the CFL conditions is violated
%     if dt*max(abs(vx_new))/dx >= 1.0 || dt*2*nu/dx^2 >= 1.0
%         fprintf('CFL condition violated, ii=%d\n',ii)
%         return
%     end
    
    % will stop running script if there are any nans in the velocity array
    nan_check = isnan(vx_new);
    
    if sum(nan_check) ~= 0
        unstable = 1;
        fprintf('unstable, ii=%d\n',ii)
        return
%     elseif rms(n_new)>=1.0e20
%         fprintf('Density value seems large.\n')
%         return
    end

    % plot loop; every 1/5 of iterations
    if plots && mod(ii,round(nmax/5))==0
        fprintf('***--------------------***\n')
        fprintf('ii=%d, count=%d\n', [ii count])
        fprintf('dt=%ds\n', dt)
        fprintf('total time=%ds\n', dt*ii)
        fprintf('simulation time %d\n', toc(timerVal))
        fprintf('Current rms tol calc %d\n', rms(vx - vx_new))
%         fprintf('total number of particles %e\n', trapz(nxax,n_new))
%         fprintf('particle balance check %e\n', trapz(nxax,n) - trapz(nxax,n_new))
%         fprintf('density source and flux balance check %e\n', rflux -...
%             trapz(vxax,source_avg*ns_mult))
%         if dt == cfl_fact*(dx^2)/(2.0*nu)
%             fprintf('Diffusive CFL condition\n')
%         elseif dt == cfl_fact*dx/max(abs(vx_new))
%             fprintf('Convective CFL condition\n')
%         end
        figure(1)
        set(gcf,'Position',[563 925 560 420])
%         semilogy(nxax(2:npts-1),n_new(2:npts-1),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        plot(nxax,n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax+ndx(1)) max(nxax-ndx(end))])
        hold on
        if MMS
            plot(nxax,ex_soln,'--','DisplayName',['exact = ' num2str(double(ii)*dt) ' s'])
        end
%         semilogy(nxax(2:npts-1),n_new_exp(2:npts-1),'--','DisplayName',['(imp)time = ' num2str(double(ii)*dt) ' s'])
        figure(2)
        set(gcf,'Position',[7 925 560 420])
        if ~MMS
            plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        elseif MMS
            plot(vxax,vx_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
            plot(vxax,ex_solu,'--','DisplayName',['exact = ' num2str(double(ii)*dt) ' s'])
        end
        xlim([min(vxax) max(vxax)])
        hold on
%         plot(vxax,vx_new_imp/cs,'--','DisplayName',['(exp)time = ' num2str(double(ii)*dt) ' s'])
        figure(3)
        set(gcf,'Position',[3 476 560 420])
        plot(vxax(2:npts-2),(vx_source(2:npts-2)),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        if MMS
            hold on
            plot(vxax,(source_stag(n_new,1,nx/2,nx/2,1,npts,ndx)),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        end
        xlim([min(vxax) max(vxax)])
        hold on
        figure(4)
        set(gcf,'Position',[563 476 560 420])
        plot(nxax(2:npts-1),n_source(2:npts-1)*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlabel('Position (m)','Fontsize',16)
        ylabel('Density source ms^{-1}','Fontsize',16)
        legend('show','Location','northwest')
        hold on
        if ~MMS
            figure(5)
            plot(vxax(2:npts-1),pf_source(2:npts-1)*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
            xlabel('Position (m)','Fontsize',16)
            ylabel('Ponderomotive source (ms^{-1})','Fontsize',16)
            legend('show','Location','northwest')
            hold on
            figure(6)
            subplot(2,1,1)
            plot(xax,real(rf_ez),'DisplayName',['Re[E_{||}], time = ' num2str(double(ii)*dt) ' s'])
            set(gca, 'XTickLabel', [])
            ylabel('E_{||} (Vm^{-1})','Fontsize',16)
            hold on
            plot(xax,imag(rf_ez),'DisplayName',['Im[E_{||}], time = ' num2str(double(ii)*dt) ' s'])
            set(gca, 'XTickLabel', [])
            legend('show','Location','northwest')
            hold on
            subplot(2,1,2)
            plot(xax,abs(rf_ez).^2,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
            xlabel('Position (m)','Fontsize',16)
            ylabel('|E_{||}|^2 (V^2m^{-2})','Fontsize',16)
            legend('show','Location','northwest')
            hold on
        end
        count = count + 1;
        
    elseif mod(ii,53)==0 || unstable
        
        transport.dt = dt;
        transport.n_source = n_source;
        transport.vx_source = vx_source;
        transport.pf_source = pf_source;
        transport.rf_ez = rf_ez;
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
        transport.xmin = xmin;
        transport.xmax = xmax;
        transport.nu = nu;
        
%         save('/Volumes/DATA/LAPD/matlab/coupled_transport.mat','-struct','transport');
        filename = strcat('C:\Users\c3149416\Documents\coupled_transport_',num2str(ii),'.mat');
        save(filename,'-struct','transport');
        
        continue
        
    end
    
end

if MMS
    l_infn = norm(ex_soln - n_new, Inf);
    l_twon = rms(ex_soln - n_new);
    l_infu = norm(ex_solu - vx_new, Inf);
    l_twou = rms(ex_solu - vx_new);
end

fprintf('***--------------------***\n')
fprintf('ii=%d, count=%d\n', [ii count])
fprintf('dt=%ds\n', dt)
fprintf('total time=%ds\n', dt*ii)
fprintf('simulation time %d\n', toc(timerVal))

        
%%

figure(1)
set(gcf,'Position',[563 925 560 420])
% semilogy(nxax(2:npts-1),n_new(2:npts-1),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
plot(nxax,n_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
hold on
if MMS
    plot(nxax,ex_soln,'--','DisplayName',['exact = ' num2str(double(ii)*dt) ' s'])
end
% semilogy(nxax(2:npts-1),n_new_exp(2:npts-1),'--','DisplayName',['(imp)time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density m^{-3}','Fontsize',16)
legend('show','Location','west')
hold off

figure(2)
set(gcf,'Position',[7 925 560 420])
hold on
if ~MMS
    plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
elseif MMS
    plot(vxax,vx_new,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
    plot(vxax,ex_solu,'--','DisplayName',['exact = ' num2str(double(ii)*dt) ' s'])
end
% plot(vxax,vx_new_imp/cs,'--','DisplayName',['(exp)time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Mach number (v/c_s)','Fontsize',16)
legend('show','Location','southeast')
hold off

figure(3)
set(gcf,'Position',[3 476 560 420])
plot(vxax(2:npts-2),vx_source(2:npts-2),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
if MMS
    hold on
    plot(vxax,(source_stag(n_new,1,nx/2,nx/2,1,npts,ndx)),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
end
xlabel('Position (m)','Fontsize',16)
ylabel('Velocity source (ms^{-1})','Fontsize',16)
legend('show','Location','northwest')
hold off

figure(4)
set(gcf,'Position',[563 476 560 420])
plot(nxax(2:npts-1),n_source(2:npts-1)*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density source (ms^{-1})','Fontsize',16)
legend('show','Location','northwest')
hold off

if ~MMS
    figure(5)
    plot(vxax(2:npts-1),pf_source(2:npts-1)*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
    xlabel('Position (m)','Fontsize',16)
    ylabel('Ponderomotive source (ms^{-1})','Fontsize',16)
    legend('show','Location','northwest')
    hold off

    figure(6)
    subplot(2,1,1)
    plot(xax,real(rf_ez),'DisplayName',['Re[E_{||}], time = ' num2str(double(ii)*dt) ' s'])
    set(gca, 'XTickLabel', [])
    ylabel('E_{||} (Vm^{-1})','Fontsize',16)
    hold on
    plot(xax,imag(rf_ez),'DisplayName',['Im[E_{||}], time = ' num2str(double(ii)*dt) ' s'])
    set(gca, 'XTickLabel', [])
    legend('show','Location','northwest')
    hold on
    subplot(2,1,2)
    plot(xax,abs(rf_ez).^2,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
    xlabel('Position (m)','Fontsize',16)
    ylabel('|E_{||}|^2 (V^2m^{-2})','Fontsize',16)
    legend('show','Location','northwest')
    hold on
end




%%

% tax = linspace(0,nmax*dt,plot_num+2);

% % for ii=1:nmax
% for jj=1:npts
%     pressure(:,jj) = (Te + Ti)*n_mat(:,jj)*e;
% end
% for jj=1:npts-1
%     pressure_av(:,jj) = 0.5*(pressure(:,jj) + pressure(:,jj+1));
%     pressure_mat(:,jj) = pressure_av(:,jj) + (1/2)*0.5*(n_mat(:,jj+1)+n_mat(:,jj))*m.*(vx_mat(:,jj).^2);
% end
% % end
% 
% figure(7)
% % levels = linspace((min(vx_mat(:))/(vx_init(2))),(max(vx_mat(:))/max(vx_init)),100);
% % levels = levels/cs;
% levels = linspace(-2.5e-2,0,100);
% % set(gca,'colorscale','log')
% contourf(vxax(2:npts-1),tax,(vx_mat(1:plot_num+2,2:npts-1) - vx_init(2:npts-1))/cs,...
%     'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar;
% 
% 
% 
% figure(8)
% % levels = linspace(round(min(n_mat(:)),-3),round(max(n_mat(:)),-3),25);
% % levels = linspace(min(n_mat(:)),max(n_mat(:)),100);
% levels = linspace(-3.0e15,3.0e15,100);
% % set(gca,'colorscale','log')
% contourf(nxax(2:npts-1),tax,n_mat(1:plot_num+2,2:npts-1) - n_init(2:npts-1),...
%     'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar

% figure(10)
% % levels = linspace((min(pressure_mat(:))),(max(pressure_mat(:))),25);
% levels = linspace(min(pressure_mat(:)),max(pressure_mat(:)),25);
% contourf(vxax,tax,pressure_mat,levels,'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar

%%

function [ans] = grad(n,dx,npts)
    ans = (n(2:npts) - n(1:npts-1))./dx;
end

function [ans] = grad2(n,dx,npts)
    cen_diff = (n(3:npts) - n(1:npts-2))/(2.0*dx);
    fwd_diff = (-3*n(1) + 4*n(2) - n(3))/(2.0*dx);
    bwd_diff = (3*n(npts) - 4*n(npts-1) + n(npts-2))/(2.0*dx);
    ans = [fwd_diff, cen_diff, bwd_diff];
end

function [ans] = avg(n,npts)
    ans = (n(2:npts) + n(1:npts-1))/2.0;
end

function [ans] = pond_source(me,mp,omega,q,Efield,dx)
    pond_const = (1.0/4.0)*((q^2)/(me*omega^2));
    for ii = 2:length(dx)
        Ediff(1,ii-1) = -(dx(1,ii)/(dx(1,ii-1)*(dx(1,ii) + dx(1,ii-1))))*Efield(1,ii-1) +...
            ((dx(1,ii) - dx(1,ii-1))/(dx(1,ii)*dx(1,ii-1)))*Efield(1,ii) +...
            (dx(1,ii-1)/(dx(1,ii)*(dx(1,ii) + dx(1,ii-1))))*Efield(1,ii+1);
    end
    ans = (1.0/mp)*pond_const*Ediff;
end

function [ans] = source_stag(n,q,Te,Ti,m,npts,dx)
    ans = -((Te + Ti)*q./(m*avg(n,npts))).*(grad(n,dx,npts));
end

function [ans] = source_col(n,q,Te,Ti,m,npts,dx)
    ans = -((Te + Ti)*q./(m*n)).*(grad2(n,dx,npts));
end

function [ans] = mms_source_cont(om,nx,knx,nxax,dt,ii,u,kux,ux,vxax,n)
    dndt = nx*om*cos(knx*nxax.^2 + om*dt*ii);
    dnudx = 2*knx*nx*nxax.*cos(knx*nxax.^2 + om*dt*ii).*u -...
        2*kux*ux.*vxax.*sin(kux*vxax.^2 + om*dt*ii).*n;
    ans = dndt + dnudx;
end

function [ans] = mms_source_mom(om,ux,kux,vxax,dt,ii,nu,u,nxax,knx,nx,n,npts)
    dudt = -om*ux*sin(kux*vxax.^2 + om*dt*ii);
    dudx = -2.0*kux*ux*vxax.*sin(kux*vxax.^2 + om*dt*ii);
    d2udx = -2.0*kux*ux*sin(kux*vxax.^2 + om*dt*ii) -...
        4.0*kux^2*ux*vxax.^2.*cos(kux*vxax.^2 + om*dt*ii);
    dndx = 2.0*knx*nx*nxax.*cos(knx*nxax.^2 + om*dt*ii);
    dndx = interp1(nxax,dndx,vxax);
    ans = dudt + u.*dudx - nu*d2udx + ((nx/2 + nx/2)./avg(n,npts)).*dndx;
end





