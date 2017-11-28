%------------------------------------------%
% ODE solve for perturbed density          %
% DVE 2015 eq 25                           %
% rlbarnett c3149416 281117                %
%------------------------------------------%

A = @(x) get_vdperp1e(x)*rot(1,1) + get_vdperp2e(x)*rot(2,1) + get_vparae(x)*rot(3,1);

B = @(x) rot(3,1)*get_gradvparaex(x) + 1i*(-om + get_vdperp1e(x)*(ky*rot(1,2) + kz*rot(1,3)) +...
    get_vdperp2e(x)*(ky*rot(2,2) + kz*rot(2,3)) + (get_vparae(x)/2.0)*(ky*rot(3,2) + kz*rot(3,3)));

get_v1ex = @(x) interp1(xax,v1e(:,1),x);
get_v1ey = @(x) interp1(xax,v1e(:,2),x);
get_v1ez = @(x) interp1(xax,v1e(:,3),x);
get_gradv1ex = @(x) gradient(get_v1x(x),dx);

C = @(x) get_N0e(x)*(get_gradv1e(x) + 1i*ky*get_v1ey(x) + 1i*kz*get_v1ez(x)) +...
    get_v1ex(x)*get_gradN0ex(x) + get_N0ex(x)*(get_v1ey(x)*lamby + get_v1ez(x)*lambz);


    