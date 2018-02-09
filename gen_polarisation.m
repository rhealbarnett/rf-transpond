%-----------------------------------------%
% Calculate possible wave polarisation    %
% rlbarnett c3149416, 100118              %
%-----------------------------------------%

%--
% initialise arrays for eigenvectors and eigenvalues
eigvec = zeros(3,npts,4);
eigval = zeros(1,npts,4);

exezkx1 = zeros(1,npts);
exezkx2 = zeros(1,npts);
exezkx3 = zeros(1,npts);
exezkx4 = zeros(1,npts);

eyezkx1 = zeros(1,npts);
eyezkx2 = zeros(1,npts);
eyezkx3 = zeros(1,npts);
eyezkx4 = zeros(1,npts);

exeykx1 = zeros(1,npts);
exeykx2 = zeros(1,npts);
exeykx3 = zeros(1,npts);
exeykx4 = zeros(1,npts);

ezeykx1 = zeros(1,npts);
ezeykx2 = zeros(1,npts);
ezeykx3 = zeros(1,npts);
ezeykx4 = zeros(1,npts);



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
%         subbed = subs(wave_eq,kperp,kperp_arr(ii,kk));

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
        exezkx1(ii) = eigvec(1,ii,1)/eigvec(3,ii,1);
        exezkx2(ii) = eigvec(1,ii,2)/eigvec(3,ii,2);
        exezkx3(ii) = eigvec(1,ii,3)/eigvec(3,ii,3);
        exezkx4(ii) = eigvec(1,ii,4)/eigvec(3,ii,4);
        
        eyezkx1(ii) = eigvec(2,ii,1)/eigvec(3,ii,1);
        eyezkx2(ii) = eigvec(2,ii,2)/eigvec(3,ii,2);
        eyezkx3(ii) = eigvec(2,ii,3)/eigvec(3,ii,3);
        eyezkx4(ii) = eigvec(2,ii,4)/eigvec(3,ii,4);
        
        exeykx1(ii) = eigvec(1,ii,1)/eigvec(2,ii,1);
        exeykx2(ii) = eigvec(1,ii,2)/eigvec(2,ii,2);
        exeykx3(ii) = eigvec(1,ii,3)/eigvec(2,ii,3);
        exeykx4(ii) = eigvec(1,ii,4)/eigvec(2,ii,4);
        
        ezeykx1(ii) = eigvec(3,ii,1)/eigvec(2,ii,1);
        ezeykx2(ii) = eigvec(3,ii,2)/eigvec(2,ii,2);
        ezeykx3(ii) = eigvec(3,ii,3)/eigvec(2,ii,3);
        ezeykx4(ii) = eigvec(3,ii,4)/eigvec(2,ii,4);       

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

%%
%--
% plot polarisation ratios
figure(20)
set(gcf,'Position',[4 68 857 886])

subplot(4,1,1)
plot(log10(N0(1:imme)),real(exezkx1(1:imme)),'.k')
hold on
plot(log10(N0(1:imme)),imag(exezkx1(1:imme)),'.r')
plot(log10(N0(1:imme)),real(exezkx3(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(exezkx3(1:imme)),'dr','MarkerSize',3)
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ex/Ez] 1', 'Im[Ex/Ez] 1','Re[Ex/Ez] 3', 'Im[Ex/Ez] 3')
set(gca, 'XTickLabel', [])
hold off

subplot(4,1,2)
plot(log10(N0(1:imme)),real(exezkx2(1:imme)),'.k')
hold on
plot(log10(N0(1:imme)),imag(exezkx2(1:imme)),'.r')
plot(log10(N0(1:imme)),real(exezkx4(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(exezkx4(1:imme)),'dr','MarkerSize',3)
ylabel('Amplitude (Vm$^{-1}$)')
xlabel('log$_{10}|$N$_0|$','Fontsize',16)
legend('Re[Ex/Ez] 2', 'Im[Ex/Ez] 2','Re[Ex/Ez] 4', 'Im[Ex/Ez] 4')
hold off

subplot(4,1,3)
plot(log10(N0(1:imme)),real(eyezkx1(1:imme)),'.k')
hold on
plot(log10(N0(1:imme)),imag(eyezkx1(1:imme)),'.r')
plot(log10(N0(1:imme)),real(eyezkx3(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(eyezkx3(1:imme)),'dr','MarkerSize',3)
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey/Ez] 1', 'Im[Ey/Ez] 1','Re[Ey/Ez] 3', 'Im[Ey/Ez] 3')
set(gca, 'XTickLabel', [])
hold off

subplot(4,1,4)
plot(log10(N0(1:imme)),real(eyezkx2(1:imme)),'.k')
hold on
plot(log10(N0(1:imme)),imag(eyezkx2(1:imme)),'.r')
plot(log10(N0(1:imme)),real(eyezkx4(1:imme)),'dk','MarkerSize',3)
plot(log10(N0(1:imme)),imag(eyezkx4(1:imme)),'dr','MarkerSize',3)
ylabel('Amplitude (Vm$^{-1}$)')
legend('Re[Ey/Ez] 2', 'Im[Ey/Ez] 2','Re[Ey/Ez] 4', 'Im[Ey/Ez] 4')
set(gca, 'XTickLabel', [])
hold off