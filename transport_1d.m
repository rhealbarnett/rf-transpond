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

