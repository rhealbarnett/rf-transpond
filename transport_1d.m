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

upwind = NaN;
central = NaN;

promptDiff = 'Differencing scheme? (upwind or central) ';
scheme = input(promptDiff, 's');

if strcmp(scheme,'upwind')
    upwind = 1;
    central = 0;
elseif strcmp(scheme, 'central')
    upwind = 0;
    central = 1;
end



%%
%--------------------------------------------------------------------------------------------------------------%
% select boundary conditions -- start with simple dirichlet
% what do you need?
% 1: to be able to specify the type of BC -- string
% 2: to be able to then set the value -- float

ldirichlet = NaN;
rdirichlet = NaN;
rnuemann = NaN;
lneumann = NaN;
periodic = NaN;

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

%%












