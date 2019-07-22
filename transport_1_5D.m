%--------------------------------------------------------------------------------------------------------------%
% 1.5D solve for transport equations (continuity and momentum)                                                                           %
% NON-CONSERVATIVE FORMS                                                                                                                                                                                %
% staggered n and vz grids                                                                                     %
% momentum eqn upwind convection, central differenced diffusion (IMEX)                                         %
% continuity eqn first order upwind (explicit, flux selecting)                                                 %
% ghost points on density 
% nabla = (kapx, kapy, dz)
% rlbarnett c3149416 190723     
%--------------------------------------------------------------------------------------------------------------%
%%

