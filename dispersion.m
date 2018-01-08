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

%--
% initialise arrays for eigenvectors and eigenvalues
evec = zeros(3,npts,4);
eval = zeros(1,npts,4);

for kk=1:4
    for ii=1:npts

        %--
        % recalculate wave equation matrix again 
        % seems super unnessesary to have this again, probably can work
        % something better than this
        wave_eq = a - we_rhs(:,:,ii);

        %--
        % substitute kx values into the matrix
        subbed = subs(wave_eq,kx,kx_arr(ii,kk));

        %--
        % calculate eigenvalues and eigenvectors
        [vecs,vals] = eig(subbed);

        %--
        % find the index of the minimum eigenvalue
        mineig = find(double(diag(vals)) == min(double(diag(vals))));

        %--
        % take eigenvector associated with the min eigenvalue
        evec(:,ii,kk) = vecs(:,mineig);
        eval(1,ii,kk) = vals(mineig,1);

    end
end

%--
% plot kx's
figure(5)
plot(xax, real(k1),'.')

hold on

plot(xax, real(k2),'.')
plot(xax, real(k3),'o')
plot(xax, real(k4),'o')
legend('k1', 'k2', 'k3', 'k4')
% ylim([-2.5, 2.5])
xlabel('Position, x ($m$)')
ylabel('Real(kx) (/m)')

hold off

figure(6)
plot(xax, imag(k1),'.')

hold on

plot(xax, imag(k2),'.')
plot(xax, imag(k3),'o')
plot(xax, imag(k4),'o')
legend('k1', 'k2', 'k3', 'k4')
% ylim([-10.0, 10.0])
xlabel('Position, x ($m$)')
ylabel('Imag(kx) (/m)')

hold off

figure(10)

subplot(3,1,1)
plot(xax,real(evec(1,:,1)),'.')
hold on
plot(xax,imag(evec(1,:,1)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex]', 'Im[Ex]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,2)
plot(xax,real(evec(2,:,1)),'.')
hold on
plot(xax,imag(evec(2,:,1)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey]', 'Im[Ey]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,3)
plot(xax,real(evec(3,:,1)),'.')
hold on
plot(xax,imag(evec(3,:,1)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ez]', 'Im[Ez]')
hold off

figure(11)

subplot(3,1,1)
plot(xax,real(evec(1,:,2)),'.')
hold on
plot(xax,imag(evec(1,:,2)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex]', 'Im[Ex]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,2)
plot(xax,real(evec(2,:,2)),'.')
hold on
plot(xax,imag(evec(2,:,2)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey]', 'Im[Ey]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,3)
plot(xax,real(evec(3,:,2)),'.')
hold on
plot(xax,imag(evec(3,:,2)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ez]', 'Im[Ez]')
hold off

figure(12)

subplot(3,1,1)
plot(xax,real(evec(1,:,3)),'.')
hold on
plot(xax,imag(evec(1,:,3)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex]', 'Im[Ex]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,2)
plot(xax,real(evec(2,:,3)),'.')
hold on
plot(xax,imag(evec(2,:,3)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey]', 'Im[Ey]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,3)
plot(xax,real(evec(3,:,3)),'.')
hold on
plot(xax,imag(evec(3,:,3)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ez]', 'Im[Ez]')
hold off

figure(13)

subplot(3,1,1)
plot(xax,real(evec(1,:,4)),'.')
hold on
plot(xax,imag(evec(1,:,4)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex]', 'Im[Ex]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,2)
plot(xax,real(evec(2,:,4)),'.')
hold on
plot(xax,imag(evec(2,:,4)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey]', 'Im[Ey]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,3)
plot(xax,real(evec(3,:,4)),'.')
hold on
plot(xax,imag(evec(3,:,4)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ez]', 'Im[Ez]')
hold off





