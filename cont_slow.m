%------------------------------------------%
% ODE solve for parallel velocity          %
% DVE 2015 eq 23                           %
% rlbarnett c3149416 151117                %
%------------------------------------------%

%------------------------------------------%
% define variables inline as anonymous     %
% functions -- avoids messy call to bvp4c  %
% function later on                        %
%------------------------------------------%

%--
% interpolate N0e over arbitrary x range
get_N0e = @(x) interp1(xax,N0e,x);

%--
% calculate A coefficient 
A = @(x) get_N0e(x).*rot(3,1);

%--
% interpolate grad*N0e over arbitrary x range
get_gradN0ex = @(x) gradient(get_N0e(x),dx);

%--
% calculate B coefficient
B = @(x) get_gradN0ex(x).*rot(3,1) + get_N0e(x).*(3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz);

%--
% interpolate N1e, v1e and perpendicular drift velocities over arbitrary x
% range
get_N1e = @(x) interp1(xax,N1e,x);
get_v1e = @(x) interp1(xax,v1e,x);
get_vdperp1e = @(x) interp1(xax,vd_perp1e,x);
get_vdperp2e = @(x) interp1(xax,vd_perp2e,x);

%--
% calculate C coefficient
C = @(x) get_gradN0ex(x).*(get_vdperp1e(x).*rot(1,1) + get_vdperp2e(x).*rot(2,1)) + get_N0e(x).*(get_vdperp1e(x).*(rot(1,2)*lamby +...
    rot(1,3)*lambz) + get_vdperp2e(x).*(rot(2,2)*lamby + rot(2,3)*lambz)) + (1.0 / 2.0).*real(gradient(conj(get_N1e(x)).*get_v1e(x),dx));

%--
% define the boundary conditions as anonymous
% functions
bound = @(ya,yb) yb - vt;

%--
% initial guess for boundary value problem solution
solinit = bvpinit(xax,vt);

%--
% call to ode_solve, inputs A(x), B(x) & C(x), bound and solinit
% outputs sol
ode_solve;

%--
% label the solutions for use in later codes
v_parae = sol.y;
gradv_paraex = sol.yp;




