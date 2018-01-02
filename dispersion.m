%-----------------------------------------%
% Find kx quartic (determinant)           %
% Solve for roots in kx                   %
% values/parameters from van eester 2015  %
% rlbarnett c3149416, 170817              %
%-----------------------------------------%

syms kx

%--
% initialise kx roots arrays, ensure they are complex
kx_arr = zeros(npts, 4);
kx_arr = complex(kx_arr);

%-- 
% k matrix
a11 = ky^2 + kz^2;
a12 = -ky*kx;
a13 = -kz*kx;
a21 = -ky*kx;
a22 = kz^2 + kx^2;
a23 = -ky*kz;
a31 = -kz*kx;
a32 = -ky*kz;
a33 = ky^2 + kx^2;

a = [[a11, a12, a13]
    [a21, a22, a23]
    [a31, a32, a33]];

%-- 
% loop through density values

for ii = 1:npts

    %--
    % wave equation rhs
    we_rhs = k0^2*cpdt;

    %--
    % set wave equation to zero to find determinant
    wave_eq = a - we_rhs(:,:,ii);

    %--
    % find kx's
    kx_quart = det(wave_eq);
    kx_coeffs = coeffs(kx_quart, 'All');
    kx_coeffs_roots = roots(kx_coeffs);
    kx_arr(ii,:) = kx_coeffs_roots;
    
end

%--
% might not be this simple -- kx root values are likely not 'ordered' (MS)
k1 = kx_arr(:,1);
k2 = kx_arr(:,2);
k3 = kx_arr(:,3);
k4 = kx_arr(:,4);

evec1 = zeros(3,npts);
evec2 = zeros(3,npts);
evec3 = zeros(3,npts);
evec4 = zeros(3,npts);

eval1 = zeros(1,npts);
eval2 = zeros(1,npts);
eval3 = zeros(1,npts);
eval4 = zeros(1,npts);

for ii=1:npts
    
    wave_eq = a - we_rhs(:,:,ii);
    subbed = subs(wave_eq,kx,k1(ii));
    [vecs,vals] = eig(subbed);
    minimum = min(double(diag(vals)));
    mineig = find(double(diag(vals)) == minimum);
    evec1(:,ii) = vecs(:,mineig);
    eval1(1,ii) = vals(mineig,1);
    
end

%--
% plot kx's
figure(5)
plot(xax, real(k1),'+')

hold on

plot(xax, real(k2),'*')
plot(xax, real(k3),'o')
plot(xax, real(k4),'^')
legend('k1', 'k2', 'k3', 'k4')
% ylim([-2.5, 2.5])
xlabel('Position, x ($m$)')
ylabel('Real(kx) (/m)')

hold off

figure(6)
plot(xax, imag(k1),'+')

hold on

plot(xax, imag(k2),'*')
plot(xax, imag(k3),'o')
plot(xax, imag(k4),'^')
legend('k1', 'k2', 'k3', 'k4')
% ylim([-10.0, 10.0])
xlabel('Position, x ($m$)')
ylabel('Imag(kx) (/m)')

hold off










