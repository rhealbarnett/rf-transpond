%------------------------------------------%
% ODE test solve                           %
% DVE 2015 eq 23                           %
% rlbarnett c3149416 151117                %
%------------------------------------------%

Av_para = rot(3,1)*N0e;
Bv_para = rot(3,1)*gradN0ex + (3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz)*N0e;
Cv_para = (vd_perp1e*rot(1,1) + vd_perp2e*rot(2,1)).*gradN0ex + (vd_perp1e*(rot(1,2)*lamby + rot(1,3)*lambz)...
    + vd_perp2e*(rot(2,2)*lamby + rot(2,3)*lambz)).*N0e + (1.0 / 2.0)*real(gradient(conj(N1e).*v1e,dx));

solinit = bvpinit(xax, [0]);

sol = bvp4c(@(x,v_para)ode(x,v_para,Av_para, Bv_para, Cv_para),@(ya,yb)bcs(ya,yb,vt),solinit);



function [ans] = get_N0e(x)

    ans = interp1(xax,N0e,x);

end

function [ans] = get_gradN0ex(x,dx)

    ans = gradient(@get_N0e,dx);

end

function [ans] = B(x)

    ans = rot(3,1)*get_gradN0ex(x) + (3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz)*get_N0e(x);

end

function dydx = ode(x,y)
    
    dydx = -(B(x)*y + C(x))/A(x);
    
end

% function dv_paradx = ode(x,v_para, Av_para, Bv_para, Cv_para)
%     
%     dv_paradx = -(Bv_para*v_para + Cv_para)/Av_para;
% %     dydx = [y(2); -abs(y(1))];
%     
% end

function res = bcs(ya,yb,vt)

    res = [ya - vt];%; yb(1) - 0.3679];
%     res = [ya(1); yb(1) + 2];
    
end




