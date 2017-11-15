%------------------------------------------%
% ODE test solve                           %
% DVE 2015 eq 23                           %
% rlbarnett c3149416 151117                %
%------------------------------------------%



Av_para = rot(3,1)*N0e;
Bv_para = rot(3,1)*gradN0e + (3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz)*N0e;
Cv_para = (vd_perp1*rot(1,1) + vd_perp2*rot(2,1))*gradN0e + (vd_perp1*(rot(1,2)*lamby + rot(1,3)*lambz)...
    + vd_perp2*(rot(2,2)*lamby + rot(2,3)*lambz))*N0e + (1.0 / 2.0)*real(ddx(conj(N1e)*v1e));



