%------------------------------------------%
% ODE solve for parallel velocity          %
% DVE 2015 eq 23                           %
% rlbarnett c3149416 151117                %
%------------------------------------------%

%------------------------------------------%
%------------------IONS--------------------%
%------------------------------------------%

vT = sqrt((T_ev*abs(e)) / me(npts/2));

%--
% interpolate N0i over arbitrary x range
get_N0i = @(x) interp1(xax,N0i,x);

%--
% calculate A coefficient 
A = @(x) get_N0i(x).*rot(3,1);

%--
% interpolate grad*N0i over arbitrary x range
get_gradN0ix = @(x) gradient(get_N0i(x),dx);

%--
% calculate B coefficient
B = @(x) get_gradN0ix(x).*rot(3,1) + get_N0i(x).*(3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz);

%--
% interpolate N1i, v1i and perpendicular drift velocities over arbitrary x
% range
get_N1i = @(x) interp1(xax,N1i,x);
get_v1ix = @(x) interp1(xax,v1i(1,:),x);
get_vdperp1i = @(x) interp1(xax,vd_perp1i,x);
get_vdperp2i = @(x) interp1(xax,vd_perp2i,x);

%--
% calculate C coefficient
C = @(x) get_gradN0ix(x).*(get_vdperp1i(x).*rot(1,1) + get_vdperp2i(x).*rot(2,1)) + get_N0i(x).*(get_vdperp1i(x).*(rot(1,2)*lamby +...
    rot(1,3)*lambz) + get_vdperp2i(x).*(rot(2,2)*lamby + rot(2,3)*lambz)) + (1.0 / 2.0).*real(gradient(conj(get_N1i(x)).*get_v1ix(x),dx));

%--
% define the boundary conditions as anonymous
% functions
bound = @(ya,yb) yb - vT;

%--
% initial guess for boundary value problem solution
solinit = bvpinit(xax,vT);

%--
% call to ode_solve, inputs A(x), B(x) & C(x), bound and solinit
% outputs sol
ode_solve;

%--
% label the solutions for use in later codes
v_parai = sol.y;
gradv_paraix = sol.yp;


%------------------------------------------%
%--------------ELECTRONS-------------------%
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
get_v1ex = @(x) interp1(xax,v1e(1,:),x);
get_vdperp1e = @(x) interp1(xax,vd_perp1e,x);
get_vdperp2e = @(x) interp1(xax,vd_perp2e,x);

%--
% calculate C coefficient
C = @(x) get_gradN0ex(x).*(get_vdperp1e(x).*rot(1,1) + get_vdperp2e(x).*rot(2,1)) + get_N0e(x).*(get_vdperp1e(x).*(rot(1,2)*lamby +...
    rot(1,3)*lambz) + get_vdperp2e(x).*(rot(2,2)*lamby + rot(2,3)*lambz)) + (1.0 / 2.0).*real(gradient(conj(get_N1e(x)).*get_v1ex(x),dx));

%--
% define the boundary conditions as anonymous
% functions
ve_bound = -N0i(end)*v_parai(end)/N0e(end);
bound = @(ya,yb) yb - ve_bound;

%--
% initial guess for boundary value problem solution
solinit = bvpinit(xax,vT);

%--
% call to ode_solve, inputs A(x), B(x) & C(x), bound and solinit
% outputs sol
ode_solve;

%--
% label the solutions for use in later codes
v_parae = sol.y;
gradv_paraex = sol.yp;




