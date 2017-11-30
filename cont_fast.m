%------------------------------------------%
% ODE solve for perturbed density          %
% DVE 2015 eq 25                           %
% rlbarnett c3149416 281117                %
%------------------------------------------%

%--
% set initial boundary locations
% opposite ends of the solution space for ions and electrons, and then
% switch
bound_oldi = N1i(1);
bound_olde = N1e(end);

%------------------------------------------%
%------------------IONS--------------------%
%------------------------------------------%

% for iter=1:nmax
%--
% A coefficient; all required quantities have been defined previously
A = @(x) get_vdperp1i(x)*rot(1,1) + get_vdperp2i(x)*rot(2,1) + get_vparai(x)*rot(3,1);

%--
% B coefficient; all required quantities have been defined previously
B = @(x) rot(3,1)*get_gradvparaix(x) + 1i*(-om + get_vdperp1i(x)*(ky*rot(1,2) + kz*rot(1,3)) +...
    get_vdperp2i(x)*(ky*rot(2,2) + kz*rot(2,3)) + (get_vparai(x)/2.0)*(ky*rot(3,2) + kz*rot(3,3)));

%--
% interpolate perturbed velocity on arbitrary grid; calculate v1 derivative
get_v1ix = @(x) interp1(xax,v1i(1,:),x);
get_v1iy = @(x) interp1(xax,v1i(2,:),x);
get_v1iz = @(x) interp1(xax,v1i(3,:),x);
get_gradv1ix = @(x) gradient(get_v1ix(x),dx);

%--
% C coefficient
C = @(x) get_N0i(x)*(get_gradv1ix(x) + 1i*ky*get_v1iy(x) + 1i*kz*get_v1iz(x)) +...
    get_v1ix(x)*get_gradN0ix(x) + get_N0i(x)*(get_v1iy(x)*lamby + get_v1iz(x)*lambz);

%--
% boundary condition
% if bound_oldi == N1i(1)
%     bound_vali = N1i(end);
% elseif bound_oldi == N1i(end)
%     bound_vali = N1i(1);
% end

bound_vali = 0.0;
bound = @(ya,yb) yb - bound_vali;
% bound_oldi = bound_vali;

solinit = bvpinit(xax,0.0);

ode_solve;

N1i = deval(sol,xax);

% end
%------------------------------------------%
%--------------ELECTRONS-------------------%
%------------------------------------------%

%--
%A coefficient; all required quantities have been defined previously
A = @(x) get_vdperp1e(x)*rot(1,1) + get_vdperp2e(x)*rot(2,1) + get_vparae(x)*rot(3,1);

%--
%B coefficient; all required quantities have been defined previously
B = @(x) rot(3,1)*get_gradvparaex(x) + 1i*(-om + get_vdperp1e(x)*(ky*rot(1,2) + kz*rot(1,3)) +...
    get_vdperp2e(x)*(ky*rot(2,2) + kz*rot(2,3)) + (get_vparae(x)/2.0)*(ky*rot(3,2) + kz*rot(3,3)));

%--
%interpolate perturbed velocity on arbitrary grid; calculate v1 derivative
get_v1ex = @(x) interp1(xax,v1e(1,:),x);
get_v1ey = @(x) interp1(xax,v1e(2,:),x);
get_v1ez = @(x) interp1(xax,v1e(3,:),x);
get_gradv1ex = @(x) gradient(get_v1ex(x),dx);

%--
%C coefficient
C = @(x) get_N0e(x)*(get_gradv1ex(x) + 1i*ky*get_v1ey(x) + 1i*kz*get_v1ez(x)) +...
    get_v1ex(x)*get_gradN0ex(x) + get_N0e(x)*(get_v1ey(x)*lamby + get_v1ez(x)*lambz);

%--
%boundary condition
bound_vale = 0.0;
bound = @(ya,yb) yb - bound_vale;
% bound_oldi = bound_vali;

solinit = bvpinit(xax,0.0);

ode_solve;

N1e = deval(sol,xax);