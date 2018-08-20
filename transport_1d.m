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

transport_params_colocated;

promptGrid = 'Staggered or collocated grid? ';
grid = input(promptGrid, 's');

if strcmp(grid,'staggered')
    staggered = 1;
    collocated = 0;
else
    staggered = 0;
    collocated = 1;
end

if staggered && size(vx_new)==size(n_new)
    fprintf("Check velocity and array density sizes for staggered grid.\n")
    return
elseif collocated && size(vx_new)~=size(n_new)
    fprintf("Check velocity and array density sizes for collocated grid.\n")
    return
end

if staggered
    
    
    for jj=2:npts-2
        if ((vx(1,jj-1)+vx(1,jj))/2)>0
            nA(jj,jj) = 1.0 - alpha*vx(1,jj);
            nA(jj,jj-1) = alpha*vx(1,jj-1);
        elseif ((vx(1,jj-1)+vx(1,jj))/2)<0
            nA(jj,jj) = 1.0 + alpha*vx(1,jj-1);
            nA(jj,jj+1) = -alpha*vx(1,jj);
        end    
    end
end

upwind = NaN;
central = NaN;

promptDiff = 'Spatial differencing scheme? (upwind or central) ';
scheme = input(promptDiff, 's');

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

ldirichlet = NaN;
rdirichlet = NaN;
rneumann = NaN;
lneumann = NaN;
periodic = NaN;

%----- need a prompt here to check whether the velocity is positive or
% negative next to the boundary, as this will determine whether left or right (or both)
% BC is required

if central

    promptlBC = 'Left BC type? (dirichlet, neumann, periodic) ';
    leftBC = input(promptlBC, 's');
    if strcmp('periodic',leftBC)

    else
        promptrBC = 'Right BC type? (dirichlet or neumann) ';
        rightBC = input(promptrBC, 's');
    end

    if strcmp('dirichlet',leftBC)
        ldirichlet = 1;               % this allows the use of ldirichlet as a logical 
        lneumann = 0;
        periodic = 0;
    elseif strcmp('neumann',leftBC)
        ldirichlet = 0;                
        lneumann = 1; 
        periodic = 0;
    end

    if strcmp('dirichlet',rightBC)
        rdirichlet = 1;
        rneumann = 0;
        periodic = 0;
    elseif strcmp('neumann',rightBC)
        rneumann = 1;
        rdirichlet = 0;
        periodic = 0;
    end

    if strcmp('periodic',leftBC)
        rdirichlet = 0;
        rneumann = 0;
        ldirichlet = 0;
        lneumann = 0;
        periodic = 1;
    end

    promptlBCval = 'Left BC value? ';
    lBC_val = input(promptlBCval);
    if strcmp('periodic',leftBC)

    else
        promptrBCval = 'Right BC value? ';
        rBC_val = input(promptrBCval);
    end

    if ldirichlet && rdirichlet
        fprintf("left and right BCs are dirichlet; lBC = %f, rBC = %f.\n", lBC_val, rBC_val)
    elseif lneumann && rdirichlet
        fprintf("left BC is neumann, lBC = %f; right BC is dirichlet, rBC = %f.\n", lBC_val, rBC_val)
    elseif rneumann && ldirichlet
        fprintf("left BC is dirichlet, lBC = %f; right BC is neumann, rBC = %f.\n", lBC_val, rBC_val)
    elseif rneumann && lneumann
        fprintf("left and right BCs are dirichlet; lBC = %f, rBC = %f.\n", lBC_val, rBC_val)
    elseif periodic
        fprintf("periodic BCs\n")
    end
    
elseif upwind
    
    if vx_new(1,2) > 0
            promptlBC = 'Left BC type? (dirichlet, neumann, periodic) ';
            leftBC = input(promptlBC, 's');
            if strcmp('dirichlet',leftBC)
                ldirichlet = 1;              
                lneumann = 0;
                periodic = 0;
            elseif strcmp('neumann',leftBC)
                ldirichlet = 0;                
                lneumann = 1; 
                periodic = 0;
            end
            promptlBCval = 'Left BC value? ';
            lBC_val = input(promptlBCval);
    elseif vx_new(1,2) < 0 
            fprintf("Left BC not required\n")
            ldirichlet = 0;
            lneumann = 0;
            periodic = 0;
    end
    
    if vx_new(1,end-1) < 0
            promptrBC = 'Right BC type? (dirichlet, neumann, periodic) ';
            rightBC = input(promptrBC, 's');
            if strcmp('dirichlet',rightBC)
                rdirichlet = 1;              
                rneumann = 0;
                periodic = 0;
            elseif strcmp('neumann',rightBC)
                rdirichlet = 0;                
                rneumann = 1; 
                periodic = 0;
            end
            promptrBCval = 'Right BC value? ';
            rBC_val = input(promptrBCval);
    elseif vx_new(1,end-1) > 0 
            fprintf("Right BC not required\n")
            rdirichlet = 0;
            rneumann = 0;
            periodic = 0;
    end
    
end
%%

if upwind
    for ii=1:nmax
        
        vx = vx_new;
        n = n_new;
        
        for jj=1:npts-1

            if jj==1
                if ldirichlet
                    vx(1,1) = lBC_val;
                    vxA(1,1) = 1.0;
                elseif lneumann
                    vx(1,1) = vx(1,2) + dx*lBC_val;
                    vxA(1,1) = 1.0;
                else
                vxA(jj,jj) = 1 + alpha*vx(1,jj);
                vxA(jj,jj+1) = -alpha*vx(1,jj);
                end
            elseif jj==npts-1
                if rdirichlet
                    vx(1,end) = rBC_val;
                    vxA(end,end) = 1.0;
                elseif rneumann
                    vx(1,end) = vx(1,end-1) + dx*rBC_val;
                    vxA(end,end) = 1.0;
                else 
                    vxA(jj,jj) = 1 - alpha*vx(1,jj);
                    vxA(jj,jj-1) = alpha*vx(1,jj);
                end
            elseif vx(1,jj)>0
                vxA(jj,jj) = 1 - alpha*vx(1,jj);
                vxA(jj,jj-1) = alpha*vx(1,jj);
            elseif vx(1,jj)<0
                vxA(jj,jj) = 1 + alpha*vx(1,jj);
                vxA(jj,jj+1) = -alpha*vx(1,jj);
            end

        end
        
        vx_new = vxA*vx';
        vx_new = vx_new';
    
    end
elseif central
    for ii=1:nmax
        
        vx = vx_new;
        n = n_new;
        
        for jj=1:npts-1

            if jj==1
                if ldirichlet
                    vx(1,1) = lBC_val;
                    vxA(1,1) = 1.0;
                elseif lneumann
                    vx(1,1) = vx(1,2) + dx*lBC_val;
                    vxA(1,1) = 1.0;
                else
                    vxA(jj,jj) = 1 + alpha*vx(1,jj);
                    vxA(jj,jj+1) = -alpha*vx(1,jj);
                end
            elseif jj==npts-1
                if rdirichlet
                    vx(1,end) = rBC_val;
                    vxA(end,end) = 1.0;
                elseif rneumann
                    vx(1,end) = vx(1,end-1) + dx*rBC_val;
                    vxA(end,end) = 1.0;
                else 
                    vxA(jj,jj) = 1 - alpha*vx(1,jj);
                    vxA(jj,jj-1) = alpha*vx(1,jj);
                end
            else
                vxA(jj,jj) = (alpha/dx)*2.0*nu;
                vxA(jj,jj-1) = -(vx(1,jj-1))*alpha/2.0 - (nu*alpha)/dx + 1/2;
                vxA(jj,jj+1) = (vx(1,jj+1))*alpha/2.0 - (nu*alpha)/dx + 1/2;
            end
        end
        
        vx_new = vxA*vx';
        vx_new = vx_new';
    end
end
        
        










