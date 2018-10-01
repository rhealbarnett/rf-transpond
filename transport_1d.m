%--------------------------------------------------------------------------------------------------------------%
% solve coupled transport equations                                                                            %
% continuity and momentum eqns                                                                                 %
% CONSERVATIVE FORMS                                                                                           %
% (partial derivatives)                                                                                        %
% dvz/dt + d(vz^2/2)/dz + (1/mn)(Te+Ti)dn/dz = 0                                                               %
% dn/dt + d(n*vz)/dz = 0                                                                                       %
% staggered n and vz grids                                                                                     %
% momentum eqn central differenced                                                                             %
% continuity eqn first order upwind (flux selecting)                                                           %
% ghost points on density for mom source term                                                                  %
%       -- first order neumann = 0                                                                             %
% rlbarnett c3149416 140818                                                                                    %
%--------------------------------------------------------------------------------------------------------------%
% set up code to accept 'switches' for left and right boundary conditions
% depending on the flux.
%
% need to use something like strcmp function to compare strings
%
% maybe also need to look into actually writing this as a function with
% whatever inputs?
%--------------------------------------------------------------------------------------------------------------%
%%
%--------------------------------------------------------------------------------------------------------------%
% --select differencing scheme
% --central should only work if there is a diffusion term
% --OR maybe set up so that the advetion term is always upwinded and the
% diffusion term is central differenced?
% --maybe this should go before the BCs, because if it is central
% differenced/has diffusion term, there will need to be two BCs


% function [n, vx] = transport_1d(grid_type,spatial_scheme,...
%   ln_bound_type,ln_bound_val,rn_bound_type,rn_bound_val,...
%   lv_bound_type,lv_bound_val,rv_bound_type,rv_bound_val,const,input_file)

% run(input_file);

% import parameter file
transport_test;

% initialise velocity and density 
vx = vx_new;
n = n_new;

staggered = NaN;
collocated = NaN;
v_ldirichlet = NaN;
v_rdirichlet = NaN;
v_rneumann = NaN;
v_lneumann = NaN;
v_periodic = NaN;
n_ldirichlet = NaN;
n_rdirichlet = NaN;
n_rneumann = NaN;
n_lneumann = NaN;
n_periodic = NaN;
explicit = NaN;
implicit = NaN;

% promptGrid = 'staggered or collocated grid? ';
grid = input(grid_type, 's');

if strcmp(grid,'staggered')
    staggered = 1;
    collocated = 0;
elseif strcmp(grid,'collocated')
    staggered = 0;
    collocated = 1;
end


if (isnan(staggered)) || (isnan(collocated))
    error("Check spelling and/or type of answer for %s.\n",'"staggered or collocated grid?"')
    return
end

if staggered
    
%     promptlGhost = 'Left (ghost) BC type? (dirichlet, neumann, periodic, linear extrap) ';
    leftGhost = input(ln_bound_type, 's');
    
    if strcmp('periodic',leftGhost)

    else
%         promptrGhost = 'Right (ghost) BC type? (dirichlet, neumann, linear extrap) ';
        rightGhost = input(rn_bound_type, 's');
    end
    
    if strcmp('linear extrap',leftGhost)
        lGhost = interp1([nxax(2), nxax(3)], [n_new(2), n_new(3)],...
            nxax(1),'linear','extrap');
        n_ldirichlet = 0;              
        n_lneumann = 0;
        n_periodic = 0;
    elseif strcmp('dirichlet',leftGhost)
        n_ldirichlet = 1;              
        n_lneumann = 0;
        n_periodic = 0;
    elseif strcmp('neumann',leftGhost)
        n_ldirichlet = 0;                
        n_lneumann = 1; 
        n_periodic = 0;
    end  
    
    if (isnan(n_ldirichlet)) || (isnan(n_lneumann)) || (isnan(n_periodic))
        error("Check spelling and/or type of answer for %s.\n",...
            '"Left (ghost) BC type? (dirichlet, neumann, periodic, linear extrap)"')
        return
    end
    
    if strcmp('linear extrap',leftGhost)
    else
%         promptlnBCval = 'Left (ghost) BC value for density? ';
        lnBC_val = input(ln_bound_val);
    end
    
    if strcmp('linear extrap',rightGhost)
        rGhost = interp1([nxax(npts-2), nxax(npts-1)], [n_new(npts-2), n_new(npts-1)],...
            nxax(npts),'linear','extrap');
        n_rdirichlet = 0;                
        n_rneumann = 0; 
        n_periodic = 0;
    elseif strcmp('dirichlet',rightGhost)
        n_rdirichlet = 1;              
        n_rneumann = 0;
        n_periodic = 0;
    elseif strcmp('neumann',rightGhost)
        n_rdirichlet = 0;                
        n_rneumann = 1; 
        n_periodic = 0;
    end
    
    if strcmp('linear extrap',rightGhost)
    else
%         promptrnBCval = 'Right (ghost) BC value for density? ';
        rnBC_val = input(rn_bound_val);
    end
    
    if (isnan(n_rdirichlet)) || (isnan(n_rneumann)) || (isnan(n_periodic))
        error("Check spelling and/or type of answer for %s.\n",...
            '"Right (ghost) BC type? (dirichlet, neumann, periodic)"')
        return
    end
    
end

if collocated
    
    if vx_new(1,2) > 0
        promptlnBC = 'Left BC type for density? (dirichlet, neumann, periodic) ';
        leftnBC = input(promptlnBC, 's');
        if strcmp('dirichlet',leftnBC)
            n_ldirichlet = 1;              
            n_lneumann = 0;
            n_periodic = 0;
        elseif strcmp('neumann',leftnBC)
            n_ldirichlet = 0;                
            n_lneumann = 1; 
            n_periodic = 0;
        end
        promptlnBCval = 'Left BC value for density? ';
        lnBC_val = input(promptlnBCval);
    elseif vx_new(1,2) < 0 
        fprintf("Left BC not required on density for the given flux direction.\n")
        n_ldirichlet = 0;
        n_lneumann = 0;
        n_periodic = 0;
    end
    if vx_new(1,end-1) < 0
        promptrnBC = 'Right BC type for density? (dirichlet, neumann, periodic) ';
        rightnBC = input(promptrnBC, 's');
        if strcmp('dirichlet',rightnBC)
            n_rdirichlet = 1;              
            n_rneumann = 0;
            n_periodic = 0;
        elseif strcmp('neumann',rightnBC)
            n_rdirichlet = 0;                
            n_rneumann = 1; 
            n_periodic = 0;
        end
        promptrnBCval = 'Right BC value for density? ';
        rnBC_val = input(promptrnBCval);
    elseif vx_new(1,end-1) > 0 
        fprintf("Right BC not required on density for the given flux direction.\n")
        n_rdirichlet = 0;
        n_rneumann = 0;
        n_periodic = 0;
    end
end

%%

upwind = NaN;
central = NaN;

% promptDiff = 'Spatial differencing scheme for momentum equation? (upwind or central) ';
scheme = input(spatial_scheme, 's');

if strcmp(scheme,'upwind')
    upwind = 1;
    central = 0;
elseif strcmp(scheme, 'central')
    upwind = 0;
    central = 1;
end

if central && nu==0
    fprintf("nu==0: central difference scheme not stable.\n")
    return
end

%%
%--------------------------------------------------------------------------------------------------------------%
% select boundary conditions -- start with simple dirichlet
% what do you need?
% 1: to be able to specify the type of BC -- string
% 2: to be able to then set the value -- float


%----- need a prompt here to check whether the velocity is positive or
% negative next to the boundary, as this will determine whether left or right (or both)
% BC is required

if central

%     promptlvBC = 'Left BC type for velocity? (dirichlet, neumann, periodic) ';
    leftvBC = input(lv_bound_type, 's');
    if strcmp('periodic',leftvBC)

    else
%         promptrvBC = 'Right BC type for velocity? (dirichlet or neumann) ';
        rightvBC = input(rv_bound_type, 's');
    end

    if strcmp('dirichlet',leftvBC)
        v_ldirichlet = 1;               % this allows the use of ldirichlet as a logical 
        v_lneumann = 0;
        v_periodic = 0;
    elseif strcmp('neumann',leftvBC)
        v_ldirichlet = 0;                
        v_lneumann = 1; 
        v_periodic = 0;
    end

    if strcmp('dirichlet',rightvBC)
        v_rdirichlet = 1;
        v_rneumann = 0;
        v_periodic = 0;
    elseif strcmp('neumann',rightvBC)
        v_rneumann = 1;
        v_rdirichlet = 0;
        v_periodic = 0;
    end

    if strcmp('periodic',leftvBC)
        v_rdirichlet = 0;
        v_rneumann = 0;
        v_ldirichlet = 0;
        v_lneumann = 0;
        v_periodic = 1;
    end

%     promptlvBCval = 'Left BC value for velocity? ';
    lvBC_val = input(lv_bound_val);
    if strcmp('periodic',leftvBC)

    else
%         promptrvBCval = 'Right BC value for velocity? ';
        rvBC_val = input(rv_bound_val);
    end

    if v_ldirichlet && v_rdirichlet
        fprintf("left and right velocity BCs are dirichlet; lBC = %f, rBC = %f.\n", lvBC_val, rvBC_val)
    elseif v_lneumann && v_rdirichlet
        fprintf("left velocity BC is neumann, lBC = %f; right velocity BC is dirichlet, rBC = %f.\n", lvBC_val, rvBC_val)
    elseif v_rneumann && v_ldirichlet
        fprintf("left velocity BC is dirichlet, lBC = %f; right velocity BC is neumann, rBC = %f.\n", lvBC_val, rvBC_val)
    elseif v_rneumann && v_lneumann
        fprintf("left and right velocity BCs are dirichlet; lBC = %f, rBC = %f.\n", lvBC_val, rvBC_val)
    elseif v_periodic
        fprintf("periodic velocity BCs\n")
    end
    
elseif upwind
    
    if vx_new(1,2) > 0
            promptlvBC = 'Left BC type for velocity? (dirichlet, neumann, periodic) ';
            leftvBC = input(promptlvBC, 's');
            if strcmp('dirichlet',leftvBC)
                v_ldirichlet = 1;              
                v_lneumann = 0;
                v_periodic = 0;
            elseif strcmp('neumann',leftvBC)
                v_ldirichlet = 0;                
                v_lneumann = 1; 
                v_periodic = 0;
            end
            promptlvBCval = 'Left BC value for velocity? ';
            lvBC_val = input(promptlvBCval);
    elseif vx_new(1,2) < 0 
            fprintf("Left BC not required\n")
            v_ldirichlet = 0;
            v_lneumann = 0;
            v_periodic = 0;
    end
    
    if vx_new(1,end-1) < 0
            promptrvBC = 'Right BC type? (dirichlet, neumann, periodic) ';
            rightvBC = input(promptrvBC, 's');
            if strcmp('dirichlet',rightvBC)
                v_rdirichlet = 1;              
                v_rneumann = 0;
                v_periodic = 0;
            elseif strcmp('neumann',rightvBC)
                v_rdirichlet = 0;                
                v_rneumann = 1; 
                v_periodic = 0;
            end
            promptrvBCval = 'Right BC value? ';
            rvBC_val = input(promptrvBCval);
    elseif vx_new(1,end-1) > 0 
            fprintf("Right BC not required\n")
            v_rdirichlet = 0;
            v_rneumann = 0;
            v_periodic = 0;
    end
    
end

%%
%--------------------------------------------------------------------------------------------------------------%
% INITALISE PLOTS; INCLUDE INITIAL CONDITIONS
%--------------------------------------------------------------------------------------------------------------%

vx_source = source_stag(n_new,e,Te,Ti,m,npts,dx);

figure(1)
set(gcf,'Position',[563 925 560 420])
semilogy(nxax(2:npts-1),n_new(2:npts-1),'DisplayName',['time = 0s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density m^{-3}','Fontsize',16)
legend('show','Location','south')
hold on

figure(2)
set(gcf,'Position',[7 925 560 420])
plot(vxax,vx_new/cs,'DisplayName',['time = 0s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Mach number','Fontsize',16)
legend('show','Location','southeast')
hold on

figure(3)
set(gcf,'Position',[3 476 560 420])
semilogy(vxax,abs(vx_source*dt),'DisplayName',['time = 0s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Velocity source ms^{-1}','Fontsize',16)
legend('show','Location','northwest')
hold on

figure(4)
set(gcf,'Position',[563 476 560 420])
plot(nxax,n_source*dt,'DisplayName',['time = 0s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density source ms^{-1}','Fontsize',16)
legend('show','Location','northwest')
hold on

%%
%--------------------------------------------------------------------------------------------------------------%
% START TIME STEPPING
%--------------------------------------------------------------------------------------------------------------%

count = 1;
timerVal = tic;

vx_rms = zeros(1,520);
n_rms = zeros(1,520);

for ii=1:520
      
    % set the vectors with the old value going into the next loop
    n = n_new;
    vx = vx_new;
    
    if staggered
        
        % fill n coefficient matrix using the averaged value of the
        % velocities on surrounding grid points
        for jj=2:npts-1
            if ((vx(1,jj-1)+vx(1,jj))/2)>0
                nA(jj,jj) = - mult*vx(1,jj);
                nA(jj,jj-1) = mult*vx(1,jj-1);
            elseif ((vx(1,jj-1)+vx(1,jj))/2)<0
                nA(jj,jj) = mult*vx(1,jj-1);
                nA(jj,jj+1) = -mult*vx(1,jj);
            end    
        end


        % build full coefficient matrices
        An_exp = nI + dt*nA;
        An_imp = nI - dt*nA;
        
        % override values in top and bottom rows to reflect neumann
        % boundary conditions for the implicit calculation
        An_imp(1,1) = 1.0; An_imp(1,2) = -1.0;
        An_imp(end,end) = 1.0; An_imp(end,end-1) = -1.0;        
        
        % calculate explicit solution
        n_new_exp = An_exp*n' + dt*n_source';
        % directly override solution vector to include neumann boundary
        % conditions for explicit method
        n_new_exp(1,1) = n_new_exp(2,1);
        n_new_exp(end,1) = n_new_exp(end-1,1);
        
        % zero old rhs values for top and bottom boundary equations for
        % implicit calculation
        n(1,1) = 0.0;
        n(1,end) = 0.0;
        % implicit calculation
        n_new_imp = An_imp\(n' + dt*n_source');
        
        % transpose solution vector
        n_new = n_new_imp;
        n_new = n_new';
        
        % only used for linearly extrapolated boundary conditions
        if strcmp('linear extrap',leftGhost)
            lGhost = interp1([nxax(2), nxax(3)], [n_new(2), n_new(3)],...
            nxax(1),'linear','extrap');
            n_new(1,1) = lGhost;
        end
        if strcmp('linear extrap',rightGhost)
            rGhost = interp1([nxax(npts-2), nxax(npts-1)], [n_new(npts-2), n_new(npts-1)],...
            nxax(npts),'linear','extrap');
            n_new(1,end) = rGhost;
        end
        

%--------------------------------------------------------------------------------------------------------------%
% COLLOCATED NOT CURRENTLY IN USE
%--------------------------------------------------------------------------------------------------------------%
      
    elseif collocated
 
        for jj=1:npts-1
            if jj==1
                if n_rdirichlet || n_rneumann
                    nA(jj,jj) = 1.0 + mult*vx(1,jj);
                    nA(jj,jj+1) = -mult*vx(1,jj+1);  
                else
                    nA(jj,jj) = 1.0 + mult*vx(1,jj);
                    nA(jj,jj+1) = -mult*vx(1,jj+1);
                end
            elseif jj==npts-1
                if n_ldirichlet || n_lneumann
                    nA(jj,jj) = 1.0 - mult*vx(1,jj);
                    nA(jj,jj-1) = mult*vx(1,jj-1); 
                else
                    nA(jj,jj) = 1.0 - mult*vx(1,jj);
                    nA(jj,jj-1) = mult*vx(1,jj-1);
                end
            elseif vx(1,jj)>0
                nA(jj,jj) = 1.0 - mult*vx(1,jj);
                nA(jj,jj-1) = mult*vx(1,jj-1);
            elseif vx(1,jj)<0
                nA(jj,jj) = 1.0 + mult*vx(1,jj);
                nA(jj,jj+1) = -mult*vx(1,jj+1);
            end
        end
        
        n_new = n_new';
        
        if n_ldirichlet
            n_new(1,1) = lnBC_val;
        elseif n_lneumann
            n_new(1,1) = n_new(1,2) + dx*lnBC_val;
        elseif n_rdirichlet
            n_new(1,end) = rnBC_val;
        elseif n_rneumann
            n_new(1,end) = n_new(1,end-1) + dx*rnBC_val;
        end
    end
    
%--------------------------------------------------------------------------------------------------------------%
% COLLOCATED (ABOVE) NOT IN USE
% BELOW IS IN USE
%--------------------------------------------------------------------------------------------------------------%
    
    % fill coefficient matrices for momentum equation, positive for v>0 and
    % negative for v<0 on the convective term; differencing of the
    % diffusion term is central and not dependent on flow direction
    for jj=2:npts-2
        if vx(1,jj)>0
            vx_pos(jj,jj) = - mult*vx(1,jj);
            vx_pos(jj,jj-1) = mult*vx(1,jj);
        elseif vx(1,jj)<0
            vx_neg(jj,jj) = mult*vx(1,jj);
            vx_neg(jj,jj+1) = -mult*vx(1,jj);
        end
        vx_diff(jj,jj) = - mult*((2.0*nu)/dx);
        vx_diff(jj,jj-1) = mult*(nu/dx);
        vx_diff(jj,jj+1) = mult*(nu/dx);
    end
         
    % construct full coefficient matrix for momentum equation
    if upwind 
        vxA = vx_pos + vx_neg;
    elseif central
        vxA = vx_pos + vx_neg + vx_diff;
    end
    
    % build full coefficient matrices
    Avx_exp = vx_I + dt*vxA;
    Avx_imp = vx_I - dt*vxA;
    % override top and bottom rows to include dirichlet boundary conditions
    % for the momentum equation (explicit and implicit methods)
    Avx_exp(1,1) = 1.0; Avx_exp(end,end) = 1.0;
    Avx_imp(1,1) = 1.0; Avx_imp(end,end) = 1.0;
    
    % calculate the source term
    if staggered
        vx_source = source_stag(n,e,Te,Ti,m,npts,dx);
    elseif collocated
        vx_source = source_col(n,e,Te,Ti,m,npts-1,dx);
    end
      
    % zero the source term at the boundaries as it is not used (dirichlet
    % boundary conditions will override the source)
    vx_source(1,1) = 0.0;
    vx_source(1,end) = 0.0;
    
    Avx_exp = vx_I + vxA*dt;
    Avx_imp = vx_I - vxA*dt;
    Avx_exp(1,1) = 1.0; Avx_exp(end,end) = 1.0;
    Avx_imp(1,1) = 1.0; Avx_imp(end,end) = 1.0;
       
    % explicit calculation
    vx_new_exp = Avx_exp*vx' + dt*vx_source';
    % implicit calculation
    vx_new_imp = Avx_imp\(vx' - dt*vx_source');
    
    % transpose solution vector
    vx_new = vx_new_imp;
    vx_new = vx_new';
    
    % reset CFL condition based on the lowest dt out of the
    % convective/diffusive CFLs
    if (cfl_fact*(dx^2)/(2.0*nu))<(cfl_fact*dx/max(abs(vx_new)))
        dt = cfl_fact*(dx^2)/(2.0*nu);
    elseif (cfl_fact*(dx^2)/(2.0*nu))>(cfl_fact*dx/max(abs(vx_new)))
        dt = cfl_fact*dx/max(abs(vx_new));
    end

    % will stop running script if either of the CFL conditions is violated
    if dt*max(abs(vx_new))/dx >= 1.0 || dt*2*nu/dx^2 >= 1.0
        fprintf('CFL condition violated, ii=%d\n',ii)
        return
    end
    
    % will stop running script if there are any nans in the velocity array
    nan_check = isnan(vx_new);
    
    if sum(nan_check) ~= 0
        fprintf('unstable, ii=%d\n',ii)
        return
    end

    % plot loop; every 1/5 of iterations
    if ii==count*round(520/5)
        fprintf('***--------------------***\n')
        fprintf('ii=%d, count=%d\n', [ii count])
        fprintf('dt=%ds\n', dt)
        fprintf('total time=%ds\n', dt*ii)
        fprintf('simulation time %d\n', toc(timerVal))
        if dt == cfl_fact*(dx^2)/(2.0*nu)
            fprintf('Diffusive CFL condition\n')
        elseif dt == cfl_fact*dx/max(abs(vx_new))
            fprintf('Convective CFL condition\n')
        end
        figure(1)
        set(gcf,'Position',[563 925 560 420])
        semilogy(nxax(2:npts-1),n_new(2:npts-1),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(nxax+dx) max(nxax-dx)])
        hold on
        semilogy(nxax(2:npts-1),n_new_exp(2:npts-1),'--','DisplayName',['(imp)time = ' num2str(double(ii)*dt) ' s'])
        figure(2)
        set(gcf,'Position',[7 925 560 420])
        plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        plot(vxax,vx_new_imp/cs,'--','DisplayName',['(exp)time = ' num2str(double(ii)*dt) ' s'])
        figure(3)
        set(gcf,'Position',[3 476 560 420])
        semilogy(vxax,abs(vx_source*dt),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlim([min(vxax) max(vxax)])
        hold on
        figure(4)
        set(gcf,'Position',[563 476 560 420])
        plot(nxax,n_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
        xlabel('Position (m)','Fontsize',16)
        ylabel('Density source ms^{-1}','Fontsize',16)
        legend('show','Location','northwest')
        hold on
        count = count + 1;
    end
    
    vx_mat(ii,:) = vx_new;
    n_mat(ii,:) = n_new;
    vx_rms(1,ii) = rms(vx_new);
    n_rms(1,ii) = rms(n_new);
    
end

        
%%

figure(1)
set(gcf,'Position',[563 925 560 420])
semilogy(nxax(2:npts-1),n_new(2:npts-1),'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
hold on
semilogy(nxax(2:npts-1),n_new_exp(2:npts-1),'--','DisplayName',['(imp)time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density m^{-3}','Fontsize',16)
legend('show','Location','south')
hold off

figure(2)
set(gcf,'Position',[7 925 560 420])
plot(vxax,vx_new/cs,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
hold on
plot(vxax,vx_new_imp/cs,'--','DisplayName',['(exp)time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Mach number','Fontsize',16)
legend('show','Location','southeast')
hold off

figure(3)
set(gcf,'Position',[3 476 560 420])
plot(vxax,vx_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Velocity source ms^{-1}','Fontsize',16)
legend('show','Location','northwest')
hold off

figure(4)
set(gcf,'Position',[563 476 560 420])
plot(nxax,n_source*dt,'DisplayName',['time = ' num2str(double(ii)*dt) ' s'])
xlabel('Position (m)','Fontsize',16)
ylabel('Density source ms^{-1}','Fontsize',16)
legend('show','Location','northwest')
hold off

%%

% for jj=1:50
%     for ii=1:npts
%         pressure(jj,ii) = (Te + Ti)*n_mat(jj,ii)*e;
%     end
%     for ii=1:npts-1
%         pressure_av(jj,ii) = 0.5*(pressure(jj,ii) + pressure(jj,ii+1));
%         pressure_mat(jj,ii) = pressure_av(jj,ii) + (1/2)*0.5*(n_mat(jj,ii+1)+n_mat(jj,ii))*m*(vx_mat(jj,ii)^2);
%     end
% end
% 
% figure(7)
% levels = linspace(round(min(vx_mat(:)),-3),round(max(vx_mat(:)),-3),25);
% contourf(vxax,tax(1:50),vx_mat(1:50,:),levels,'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar
% 
% figure(8)
% % levels = linspace(round(min(n_mat(:)),-3),round(max(n_mat(:)),-3),25);
% levels = linspace(min(n_mat(:)),max(n_mat(:)),25);
% contourf(nxax(2:npts-1),tax(1,1:50),n_mat(1:50,2:npts-1),levels,'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar
% 
% figure(10)
% % levels = linspace((min(pressure_mat(:))),(max(pressure_mat(:))),25);
% levels = linspace(min(pressure_mat(:)),max(pressure_mat(:)),25);
% contourf(vxax,tax(1,1:50),pressure_mat(1:50,:),levels,'LineColor','none')
% xlabel('Position (m)','Fontsize',16); ylabel('Time (s)','Fontsize',16)
% colorbar

%%

function [ans] = grad(n,dx,npts)
    ans = (n(2:npts) - n(1:npts-1))/dx;
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

% function [ans] = pond_source(m,omega,q,Efield,dx,npts)
%     pond_const = (1.0/4.0)*((q^2)/(m*omega^2));
%     ans = (1.0/m)*pond_const*grad(Efield,dx,npts);
% end

function [ans] = source_stag(n,q,Te,Ti,m,npts,dx)
    ans = -((Te + Ti)*q./(m*avg(n,npts))).*(grad(n,dx,npts));
end

function [ans] = source_col(n,q,Te,Ti,m,npts,dx)
    ans = -((Te + Ti)*q./(m*n)).*(grad2(n,dx,npts));
end

end




