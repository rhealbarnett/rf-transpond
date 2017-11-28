%------------------------------------------%
% ODE solve for perturbed density          %
% DVE 2015 eq 25                           %
% rlbarnett c3149416 281117                %
%------------------------------------------%

A = @(x) get_vdperp1e(x)*rot(1,1) + get_vdperp2e(x)*rot(2,1) + get_vparae(x)*rot(3,1);

B = @(x) rot(3,1)*get_gradvparaex(x) + 1i*(-om + get_vdperp1e(x)*(ky*rot(1,2) + kz*rot(1,3)) +...
    get_vdperp2e(x)*(ky*rot(2,2) + kz*rot(2,3)) + (get_vparae(x)/2.0)*(ky*rot(3,2) + kz*rot(3,3)));


    