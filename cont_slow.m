%------------------------------------------%
% ODE test solve                           %
% DVE 2015 eq 23                           %
% rlbarnett c3149416 151117                %
%------------------------------------------%

% Av_para = rot(3,1)*N0e;
% Bv_para = rot(3,1)*gradN0ex + (3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz)*N0e;
% Cv_para = (vd_perp1e*rot(1,1) + vd_perp2e*rot(2,1)).*gradN0ex + (vd_perp1e*(rot(1,2)*lamby + rot(1,3)*lambz)...
%     + vd_perp2e*(rot(2,2)*lamby + rot(2,3)*lambz)).*N0e + (1.0 / 2.0)*real(gradient(conj(N1e).*v1e,dx));

% data.xax = xax;
% data.N0e = N0e;
% data.rot = rot;

get_N0e = @(x) interp1(xax,N0e,x);
A = @(x) get_N0e(x)*rot(3,1);
get_gradN0ex = @(x) gradient(get_N0e(x),dx);
B = @(x) get_gradN0ex(x)*rot(3,1) + get_N0e(x)*(3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz);
get_N1e = @(x) interp1(xax,N1e,x);
get_v1e = @(x) interp1(xax,v1e,x);
get_vdperp1e = @(x) interp1(xax,vd_perp1e,x);
get_vdperp2e = @(x) interp1(xax,vd_perp2e,x);
C = @(x) get_gradN0ex(x)*(get_vdperp1e(x)*rot(1,1) + get_vdperp2e(x)*rot(2,1)) + get_N0e(x)*(get_vdperp1e(x)*(rot(1,2)*lamby + rot(1,3)*lambz)...
        + get_vdperp2e(x)*(rot(2,2)*lamby + rot(2,3)*lambz)) + (1.0 / 2.0)*real(gradient(conj(get_N1e(x))*get_v1e(x),dx));

dydx = @(x,y) -(B(x)*y + C(x))/A(x);

solinit = bvpinit(xax, [0]);

sol = bvp4c(dydx,@(ya,yb)bcs(ya,yb,vt),solinit);
    
% function [ans] = C(x)
% 
%     ans = (VD_perp1e(x)*rot(1,1) + VD_perp2e(x)*rot(2,1)).*get_gradN0ex(x) + (VD_perp1(x)*(rot(1,2)*lamby + rot(1,3)*lambz)...
%         + VD_perp2(x)*(rot(2,2)*lamby + rot(2,3)*lambz)).*get_N0e(x) + (1.0 / 2.0)*real(gradient(conj(get_N1e(x)).*get_v1e(x),dx));
%     
% end

% function [VD_perp1e,VD_perp2e] = get_vdperp(x)
% 
%         VD_perp1e = interp1(xax,vd_perp1e,x);
%         VD_perp2e = interp1(xax,vd_perp2e,x);
% 
% end

% function [ans] = get_N1e(x)
% 
%     ans = interp1(xax,N1e,x);
% 
% end
% 
% function [ans] = get_v1e(x)
% 
%     ans = interp1(xax,v1e,x);
% 
% end

% function [ans] = B(x,rot)
% 
%     ans = rot(3,1)*get_gradN0ex(x) + (3.0 / 2.0)*(rot(3,2)*lamby + rot(3,3)*lambz)*get_N0e(x);
%     
% end
% 
% function [ans] = get_gradN0ex(x,dx)
% 
%     ans = gradient(@get_N0e,dx);
%     
% end


% function [ans] = A(x,data)
% 
%     ans = rot(3,1)*get_N0e(x,data);
%     
% end
% 
% 
% function [ans] = get_N0e(x,data)
% 
%     ans = interp1(data.xax,data.N0e,x);
%     
% end
% 
% function dydx = ode(func)
%     
%     dydx = func;
%     
% end


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




