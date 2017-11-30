%------------------------------------------%
% ODE solve for lnN0                       %
% DVE 2015 eq 24; second set of coeffs     %
% rlbarnett c3149416 211117                %
%------------------------------------------%

%------------------------------------------%
%------------------IONS--------------------%
%------------------------------------------%

%--
% interpolate parallel velocity over arbitrary grid x
get_vparai = @(x) interp1(xax,v_parai,x);

%--
% calculate square of parallel velocity
v_paraisq = v_parai.^2;

%--
% interpolate over arbitrary grid x
get_vparaisq = @(x) interp1(xax,v_paraisq,x);

%--
% A coefficient
A = @(x) (get_vparaisq(x) - vt^2)*rot(3,1) + get_vparai(x).*(get_vdperp1i(x)*rot(1,1) + get_vdperp2i(x)*rot(2,1));

%--
% B coefficient
B = @(x) (get_vparaisq(x) - vt^2)*(rot(3,2)*lamby + rot(3,3)*lambz) + get_vparai(x).*(get_vdperp1i(x)*(rot(1,2)*lamby +...
    rot(1,3)*lambz) + get_vdperp2i(x)*(rot(2,2)*lamby + rot(2,3)*lambz));

%--
% interpolate parallel velocity derivative, static potential, static Ex and a_pondx 
get_gradvparaix = @(x) interp1(xax,gradv_paraix,x);
get_staticpot = @(x) interp1(xax,static_pot,x);
get_staticex = @(x) interp1(xax,static_ex,x);
get_staticey = @(x) interp1(xax,static_ey,x);
get_staticez = @(x) interp1(xax,static_ez,x);
get_apondix = @(x) interp1(xax,a_pondix,x);
get_apondiy = @(x) interp1(xax,a_pondiy,x);
get_apondiz = @(x) interp1(xax,a_pondiz,x);

%--
% C coefficient
C = @(x) -get_vdperp1i(x).*(rot(1,1).*get_gradvparaix(x) + (get_vparai(x)./2.0)*(rot(1,2)*lamby + rot(1,3)*lambz)) -...
    get_vdperp2i(x).*(rot(2,1).*get_gradvparaix(x) + (get_vparai(x)./2.0)*(rot(2,2)*lamby + rot(2,3)*lambz)) +...
    ((e)/me)*(rot(3,1).*get_staticex(x) + rot(3,2).*get_staticey(x) + rot(3,3).*get_staticez(x)) + rot(3,1).*get_apondix(x) +...
    (get_vparai(x)./(2.0.*get_N0i(x))).*real(gradient(conj(get_N1i(x)).*get_v1ix(x),dx));

%--
% bounday condition 
bound = @(ya,yb) yb - log(Nmax);

%--
% initial guess for boundary value problem solution
test = linspace(xmin,xmax,10*npts);
solinit = bvpinit(test,log(Nmax)/2);

%--
% call to ode_solve, inputs A(x), B(x) & C(x), bound and solinit
% outputs sol
ode_solve;

%--
% label solutions
N0ilog = deval(sol,xax);
grad_N0ixlog = sol.yp;
N0i = exp(N0ilog);
grad_N0ix = gradient(N0i,dx);

%------------------------------------------%
%--------------ELECTRONS-------------------%
%------------------------------------------%

%--
% interpolate parallel velocity over arbitrary grid x
get_vparae = @(x) interp1(xax,v_parae,x);

%--
% calculate square of parallel velocity
v_paraesq = v_parae.^2;

%--
% interpolate over arbitrary grid x
get_vparaesq = @(x) interp1(xax,v_paraesq,x);

%--
% A coefficient
A = @(x) (get_vparaesq(x) - vt^2)*rot(3,1) + get_vparae(x).*(get_vdperp1e(x)*rot(1,1) + get_vdperp2e(x)*rot(2,1));

%--
% B coefficient
B = @(x) (get_vparaesq(x) - vt^2)*(rot(3,2)*lamby + rot(3,3)*lambz) + get_vparae(x).*(get_vdperp1e(x)*(rot(1,2)*lamby +...
    rot(1,3)*lambz) + get_vdperp2e(x)*(rot(2,2)*lamby + rot(2,3)*lambz));

%--
% interpolate parallel velocity derivative and a_pondx 
get_gradvparaex = @(x) interp1(xax,gradv_paraex,x);
get_apondex = @(x) interp1(xax,a_pondex,x);
get_apondey = @(x) interp1(xax,a_pondey,x);
get_apondez = @(x) interp1(xax,a_pondez,x);

%--
% C coefficient
C = @(x) -get_vdperp1e(x).*(rot(1,1).*get_gradvparaex(x) + (get_vparae(x)./2.0)*(rot(1,2)*lamby + rot(1,3)*lambz)) -...
    get_vdperp2e(x).*(rot(2,1).*get_gradvparaex(x) + (get_vparae(x)./2.0)*(rot(2,2)*lamby + rot(2,3)*lambz)) +...
    ((e)/me)*(rot(3,1).*get_staticex(x) + rot(3,2).*get_staticey(x) + rot(3,3).*get_staticez(x)) + rot(3,1).*get_apondex(x) +...
    (get_vparae(x)./(2.0.*get_N0e(x))).*real(gradient(conj(get_N1e(x)).*get_v1ex(x),dx));

%--
% bounday condition 
bound = @(ya,yb) yb - log(Nmax);

%--
% initial guess for boundary value problem solution
test = linspace(xmin,xmax,10*npts);
solinit = bvpinit(test,log(Nmax)/2);

%--
% call to ode_solve, inputs A(x), B(x) & C(x), bound and solinit
% outputs sol
ode_solve;

%--
% label solutions
N0elog = deval(sol,xax);
grad_N0exlog = sol.yp;
N0e = exp(N0elog);
grad_N0ex = gradient(N0e,dx);




