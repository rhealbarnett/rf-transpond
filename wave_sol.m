%--------------------------------%
% time independent wave solver   %
% V X V X E = k0^2.K.E           %
% rlbarnett c3149416 210917      %
%--------------------------------%

cpdt = zeros(3,3,npts);

%--
% cold plasma dielectric tensor elements
for nn=1:npts
    s = 1.0 - om_pe(nn)^2/(om^2 - om_ce^2) - om_pd(nn)^2/(om^2 - om_cd^2) - om_ph(nn)^2/(om^2 - om_ch^2);
    d = om_ce*om_pe(nn)^2/(om*(om^2 - om_ce^2)) + om_cd*om_pd(nn)^2/(om*(om^2 - om_cd^2)) + om_ch*om_ph(nn)^2/(om*(om^2 - om_ch^2));
    p = 1.0 - om_pe(nn)^2/om^2 - om_pd(nn)^2/om^2 - om_ph(nn)^2/om^2;

    %--
    % cold plasma delectric tensor
    cpdt(1,1,nn) = s;
    cpdt(2,2,nn) = s;
    cpdt(1,2,nn) = -1i*d;
    cpdt(2,1,nn) = 1i*d;
    cpdt(3,3,nn) = p;
    
    cpdt(:,:,nn) = rot'*cpdt(:,:,nn)*rot;
end
    
%--
% set up rhs matrix (multiples E field)

% syms ky kz k0 dx
% cpdt = sym('R%d%d',[3,3]);

A = zeros(3*npts, 3*npts);
% waveeq_mat = sym('O%d%d', [3*npts,3*npts]);
ii = 4;
kk = 1;

for eq1=4:3:3*(npts-1)
    
    eq2 = eq1+1;
    eq3 = eq2+1;
    
    iiex = ii;
    iiey = ii+1;
    iiez = ii+2;
    iiexm = iiex - 3;
    iieym = iiey - 3;
    iiezm = iiez - 3;
    iiexp = iiex + 3;
    iieyp = iiey + 3;
    iiezp = iiez + 3;

    %--
%     % set up periodic boundary conditions
%     if ((iiexm) & (iieym) & (iiezm)) <= 0
%         iiexm = 3*npts + iiexm;
%         iieym = 3*npts + iieym;
%         iiezm = 3*npts + iiezm;
%     end
    
    %--
%     % this loop doesn't work if it's set up like the previous one???
%     if ((iiexp) & (iieyp) & (iiezp)) > (3*npts)
%         iiexp = iiexp - 3*npts;
%         iieyp = iieyp - 3*npts;
%         iiezp = iiezp - 3*npts;
%     end          
    
    %--
    % fill matrix
    A(eq1,iiexm) = 0.0;
    A(eq1,iieym) = -1i*ky;
    A(eq1,iiezm) = -1i*kz;
    A(eq1,iiex) = 2.0*dx*(ky^2 + kz^2 - k0^2*cpdt(1,1,kk));
    A(eq1,iiey) = -2.0*dx*k0^2*cpdt(1,2,kk);
    A(eq1,iiez) = -2.0*dx*k0^2*cpdt(1,3,kk);
    A(eq1,iiexp) = 0.0;
    A(eq1,iieyp) = 1i*ky;
    A(eq1,iiezp) = 1i*kz;
    
    A(eq2,iiexm) = -1i*ky*(dx/2.0);
    A(eq2,iieym) = -1.0;
    A(eq2,iiezm) = 0.0;
    A(eq2,iiex) = -dx^2*k0^2*cpdt(2,1,kk);
    A(eq2,iiey) = dx^2*(kz^2 - k0^2*cpdt(2,2,kk)) + 2.0;
    A(eq2,iiez) = -dx^2*(ky*kz + k0^2*cpdt(2,3,kk));
    A(eq2,iiexp) = 1i*ky*(dx/2.0);
    A(eq2,iieyp) = -1.0;
    A(eq2,iiezp) = 0.0;
    
    A(eq3,iiexm) = -1i*kz*(dx/2.0);
    A(eq3,iieym) = 0.0;
    A(eq3,iiezm) = -1.0;
    A(eq3,iiex) = -dx^2*k0^2*cpdt(3,1,kk);
    A(eq3,iiey) = -dx^2*(ky*kz + k0^2*cpdt(3,2,kk));
    A(eq3,iiez) = dx^2*(ky^2 - k0^2*cpdt(3,3,kk)) + 2.0;
    A(eq3,iiexp) = 1i*kz*(dx/2.0);
    A(eq3,iieyp) = 0.0;
    A(eq3,iiezp) = -1.0;
    
    ii = ii + 3;
    kk = kk + 1;

end

%--
% metallic wall BC
A(1,1:6) = 0.0;
A(2,1:6) = 0.0;
A(3,1:6) = 0.0;
A(1,1) = 1.0;
A(2,2) = 1.0;
A(3,3) = 1.0;

A(3*npts-2,1:3) = 0.0;
A(3*npts-1,1:3) = 0.0;
A(3*npts,1:3) = 0.0;
A(3*npts-2,3*npts-2) = 1.0;
A(3*npts-1,3*npts-1) = 1.0;
A(3*npts-0,3*npts-0) = 1.0;

% waveeq_mat = sparse(waveeq_mat);

%--
% set up rhs vector
Jy = 1.0;
Jz = 1.0;
xloc = find(xax<=0.19);
xloc = xloc*3.0;
rhs = zeros(3*npts,1);
rhs(xloc(end)) = 0.0;
rhs(xloc(end)+1) = 1i*om*mu0*Jy;
rhs(xloc(end)+2) = 1i*om*mu0*Jz;

% --
% calculation solution as waveeq_mat^-1*rhs
% -- COMMENT OUT IF DOING SYMBOLIC MATRIX --%
rf_e = (A)\rhs;

rf_ex = rf_e(1:3:3*npts);
rf_ey = rf_e(2:3:3*npts);
rf_ez = rf_e(3:3:3*npts);

% ----------------------plots----------------------- %
figure(1)

subplot(3,1,1)
plot(xax, real(rf_ex), 'k')
ylabel('Amplitude (?)')

hold on

plot(xax, imag(rf_ex), '--')
set(gca, 'XTickLabel', [])
legend('Re[Ex]', 'Im[Ex]', 'Location', 'northwest')
%vline(0.19, '--r', 'Antenna')
    
hold off

subplot(3,1,2)
plot(xax, real(rf_ey), 'k')
ylabel('Amplitude (?)')

hold on

plot(xax, imag(rf_ey), '--')
set(gca, 'XTickLabel', [])
legend('Re[Ey]', 'Im[Ey]', 'Location', 'northwest')
%vline(0.19, '--r', 'Antenna')

hold off

subplot(3,1,3)
plot(xax, real(rf_ez), 'k')
ylabel('Amplitude (?)')

hold on

plot(xax, imag(rf_ez), '--')
xlabel('Position')
legend('Re[Ez]', 'Im[Ez]', 'Location', 'northwest')
%vline(0.19, '--r', 'Antenna')

hold off




























