%-----------------------------------------%
% Calculate possible wave polarisation    %
% rlbarnett c3149416, 100118              %
%-----------------------------------------%

%--
% initialise arrays for eigenvectors and eigenvalues
eigvec = zeros(3,npts,4);
eigval = zeros(1,npts,4);
exezkp11 = zeros(1,npts);
exezkp12 = zeros(1,npts);
exezkp21 = zeros(1,npts);
exezkp22 = zeros(1,npts);

for kk=1:4
    for ii=1:npts

        %--
        % recalculate wave equation matrix again 
        % seems super unnessesary to have this again, probably can work
        % something better than this
        wave_eq = a - we_rhs(:,:,ii);

        %--
        % substitute kx values into the matrix
%         subbed = subs(wave_eq,kx,kx_arr(ii,kk));
        subbed = subs(wave_eq,kperp,sf_arr(ii,kk));

        %--
        % calculate eigenvalues and eigenvectors
        subbed = eval(subbed);
        [vecs,vals] = eig(subbed);

        %--
        % find the index of the minimum eigenvalue
        mineig = find((diag(vals)) == min((diag(vals))));

        %--
        % take eigenvector associated with the min eigenvalue
        eigvec(:,ii,kk) = vecs(:,mineig);
        eigval(1,ii,kk) = vals(mineig,1);
        
        %--
        % calculate polarisation ratios
        exezkp11(ii) = eigvec(1,ii,1)/eigvec(3,ii,1);
        exezkp12(ii) = eigvec(1,ii,2)/eigvec(3,ii,2);
        exezkp21(ii) = eigvec(1,ii,3)/eigvec(3,ii,3);
        exezkp22(ii) = eigvec(1,ii,4)/eigvec(3,ii,4);

    end
end

%%
%--
% plot eigenvectors (ie polarisation)
figure(10)
set(gcf,'Position',[0   536   824   419])
suptitle('kp1')

subplot(3,1,1)
plot(log10(N0),real(eigvec(1,:,1)),'.')
hold on
plot(log10(N0),imag(eigvec(1,:,1)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex]', 'Im[Ex]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,2)
plot(log10(N0),real(eigvec(2,:,1)),'.')
hold on
plot(log10(N0),imag(eigvec(2,:,1)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey]', 'Im[Ey]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,3)
plot(log10(N0),real(eigvec(3,:,1)),'.')
hold on
plot(log10(N0),imag(eigvec(3,:,1)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
legend('Re[Ez]', 'Im[Ez]')
hold off

figure(11)
set(gcf,'Position',[859   536   824   419])
suptitle('kp1')

subplot(3,1,1)
plot(log10(N0),real(eigvec(1,:,2)),'.')
hold on
plot(log10(N0),imag(eigvec(1,:,2)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex]', 'Im[Ex]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,2)
plot(log10(N0),real(eigvec(2,:,2)),'.')
hold on
plot(log10(N0),imag(eigvec(2,:,2)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey]', 'Im[Ey]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,3)
plot(log10(N0),real(eigvec(3,:,2)),'.')
hold on
plot(log10(N0),imag(eigvec(3,:,2)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
legend('Re[Ez]', 'Im[Ez]')
hold off

figure(12)
set(gcf,'Position',[6    60   824   419])
suptitle('kp2')

subplot(3,1,1)
plot(log10(N0),real(eigvec(1,:,3)),'.')
hold on
plot(log10(N0),imag(eigvec(1,:,3)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex]', 'Im[Ex]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,2)
plot(log10(N0),real(eigvec(2,:,3)),'.')
hold on
plot(log10(N0),imag(eigvec(2,:,3)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey]', 'Im[Ey]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,3)
plot(log10(N0),real(eigvec(3,:,3)),'.')
hold on
plot(log10(N0),imag(eigvec(3,:,3)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
legend('Re[Ez]', 'Im[Ez]')
hold off

figure(13)
set(gcf,'Position',[841 71 835 420])
suptitle('kp2')

subplot(3,1,1)
plot(log10(N0),real(eigvec(1,:,4)),'.')
hold on
plot(log10(N0),imag(eigvec(1,:,4)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex]', 'Im[Ex]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,2)
plot(log10(N0),real(eigvec(2,:,4)),'.')
hold on
plot(log10(N0),imag(eigvec(2,:,4)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey]', 'Im[Ey]')
set(gca, 'XTickLabel', [])
hold off

subplot(3,1,3)
plot(log10(N0),real(eigvec(3,:,4)),'.')
hold on
plot(log10(N0),imag(eigvec(3,:,4)),'.')
ylabel('Amplitude (Vm$^{-1}$)')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
legend('Re[Ez]', 'Im[Ez]')
hold off
