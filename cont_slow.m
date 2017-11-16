%------------------------------------------%
% ODE test solve                           %
% DVE 2015 eq 23                           %
% rlbarnett c3149416 151117                %
%------------------------------------------%

Av_para = rot(3,1)*N0e;
Bv_para = rot(3,1)*gradN0ex + (3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz)*N0e;
Cv_para = (vd_perp1e*rot(1,1) + vd_perp2e*rot(2,1)).*gradN0ex + (vd_perp1e*(rot(1,2)*lamby + rot(1,3)*lambz)...
    + vd_perp2e*(rot(2,2)*lamby + rot(2,3)*lambz)).*N0e + (1.0 / 2.0)*real(gradient(conj(N1e).*v1e,dx));

bound = vt;

[x,v_para] = ode45(@(x,v_para) myode(x, v_para, xax, Av_para, Bv_para, Cv_para), xax, bound);

function dv_paradx = myode(x, v_para, xax, Av_para, Bv_para, Cv_para)

    Av_para = interp1(xax,Av_para,x);
    Bv_para = interp1(xax,Bv_para,x);
    Cv_para = interp1(xax,Cv_para,x);

    dv_paradx = -(Bv_para*v_para + Cv_para)./Av_para;

end




