%-----------------------------------------%
% Find kz quartic (determinant)           %
% Solve for roots in kz                   %
% can solve when B0 is not aligned with z %
% rlbarnett c3149416, 010218              %
%-----------------------------------------%

% ky = 0;
% kx = 20i;
syms kz

%--
% initialise kx roots arrays, ensure they are complex
kz_arr = zeros(npts, 4);
kz_arr = complex(kz_arr);

kp11 = zeros(npts,length(n_refrac));
kp12 = zeros(npts,length(n_refrac));
kp21 = zeros(npts,length(n_refrac));
kp22 = zeros(npts,length(n_refrac));
% k = kx;
% ky = 0.0;

for jj=1:length(n_refrac)
%     jj
    %-- 
    % k matrix
    a11 = ky(1,jj)^2 + kz^2;
    a12 = -ky(1,jj)*kx;
    a13 = -kz*kx;
    a21 = -ky(1,jj)*kx;
    a22 = kz^2 + kx^2;
    a23 = -ky(1,jj)*kz;
    a31 = -kz*kx;
    a32 = -ky(1,jj)*kz;
    a33 = ky(1,jj)^2 + kx^2;

    a = [[a11, a12, a13]
        [a21, a22, a23]
        [a31, a32, a33]];

    % kx_quart_arr = sym('K',[npts,1]);
    check = zeros(3, 3);

    %--
    % wave equation rhs
    we_rhs = k0^2*cpdt;

    %-- 
    % loop through density values

    for ii = 1:npts

        %--
        % set wave equation to zero to find determinant
        wave_eq = a - we_rhs(:,:,ii);

        %--
        % the determinant of the above equation will be a quartic in kx -- the
        % dispersion relation
        kz_quart = det(wave_eq);
    %     kx_quart_arr(ii,1) = kx_quart;

        %--
        % coeffs + 'All' finds the polynomial coeffients on the highest to
        % lowest order terms (ie for ax^4 + bx^3 ... etc they are ordered [a,
        % b, c, d, e]
        kz_coeffs = coeffs(kz_quart, 'All');

        %--
        % the roots function uses the polynomial coefficients, in order highest
        % to lowest, to determine the polynomial roots. 
        kz_coeffs_roots = roots(kz_coeffs);

        %--
        % store the four roots
        kz_arr(ii,:) = kz_coeffs_roots;

%         for kk=1:4
%             check(ii,kk) = vpa(subs(wave_eq,kz,kz_arr(ii,kk)));
%         end

    end

    %--
    % might not be this simple -- kx root values are likely not 'ordered' (MS)
    kp11(:,jj) = kz_arr(:,1);
    kp12(:,jj) = kz_arr(:,2);
    kp21(:,jj) = kz_arr(:,3);
    kp22(:,jj) = kz_arr(:,4);
end
%%

kpara11 = k1;
kpara12 = k3;
kpara21 = k2;
kpara22 = k4;

width = 1000;
height = 800;
x0 = 0;
y0 = 0;

indn_np = find(n_new<=1.0e17,1,'first');
indky_np = find(ky<=30,1,'last');
indn_p = find(n_new<=1.0e17,1,'first');
indky_p = find(ky==0,1,'last');
indn_b = find(n_new<=1.0e17,1,'last');
indky_b = find(ky<=20,1,'last');
indn_tm = find(n_new<=1.0e16,1,'last');
indky_tm = find(ky==0,1,'last');
indn_bn = find(n_new<=1.0e18,1,'first');
indky_bn = find(ky>=20,1,'first');

figure(1)
set(gcf,'Position',[x0 y0 width height],'Color','w')
ax1 = subplot(2,2,1);
levelsrkp11 = linspace(0,30,50);
contourf(log10(n_new),abs(ky),real(kpara11)',levelsrkp11,'Linecolor','none')
colormap(magma)
rc11=colorbar;
rc11.Ruler.TickLabelFormat = '%1.1f';
rc11.Ruler.Exponent = 1;
hold on
contour(log10(n_new),abs(ky),real(kpara11)',[0.01665,0.01665],'r','Linewidth',2);
plot(log10(n_new(indn_np)),ky(indky_np),'sr','MarkerSize',10,'Linewidth',2,...
    'MarkerFaceColor','r')
plot(log10(n_new(indn_p)),ky(indky_p),'.r','MarkerSize',40)
plot(log10(n_new(indn_b)),ky(indky_b),'*r','MarkerSize',10,'Linewidth',2)
plot(log10(n_new(indn_tm)),ky(indky_tm),'xr','MarkerSize',15,'Linewidth',3)
plot(log10(n_new(indn_bn)),ky(indky_bn),'dr','MarkerSize',10,'MarkerFaceColor','r',...
    'Linewidth',3)
set(gca,'xtick',[])
title('Re[k_{|| 11}] (m^{-1})')
ylabel('{\itk_{y}} (m^{-1})')
set(gca,'Fontsize',30,'FontName','CMU Serif')
set(gca,'Position',[0.15 0.55 0.3 0.35])
h1 = get(rc11,'Position');
text(0.02,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
    'color','white')
hold off

ax2 = subplot(2,2,2);
levelsikp11 = linspace(0,2.5,50);
contourf(log10(n_new),abs(ky),imag(kpara11)',levelsikp11,'Linecolor','none')
ic11=colorbar();
colormap(magma)
ic11.Ruler.TickLabelFormat = '%1.1f';
hold on
contour(log10(n_new),abs(ky),imag(kpara11)',[0.0177,0.0177],'r','Linewidth',2);
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Im[k_{|| 11}] (m^{-1})')
ylabel('k_{y} (m^{-1})')
set(gca,'Fontsize',30,'FontName','CMU Serif')
set(gca,'Position',[0.55 0.55 0.3 0.35])
set(gca,'colorscale','log')
h2 = get(ic11,'Position');
text(0.02,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
    'color','white')
hold off

ax3 = subplot(2,2,3);
levelsrkp21 = linspace(0,0.25,50);
contourf(log10(n_new),abs(ky),real(kpara21)',levelsrkp21,'Linecolor','none')
hold on
contourf(log10(n_new),abs(ky),real(kpara21)',[0 0],'Linecolor','k')
rc21=colorbar();
rc21.Ruler.Exponent = -1;
caxis([0 1])
rc21.Ruler.TickLabelFormat = '%1.1f';
colormap(magma)
hold on
contour(log10(n_new),abs(ky),real(kpara21)',[0.0027,0.0027],'-.w','Linewidth',3);
title('Re[k_{|| 21}] (m^{-1})')
ylabel('{\itk_{y}} (m^{-1})'),
xlabel('log_{10}n')
text(0.02,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
    'color','white')
set(gca,'Fontsize',30,'FontName','CMU Serif')
set(gca,'colorscale','log')
h = get(gca,'colorbar');
set(rc21,'Position',[h1(1) 0.12 h1(3) h1(4)])
set(gca,'Position',[0.15 0.12 0.3 0.35])

ax4 = subplot(2,2,4);
levelsikp21 = linspace(0,35,50);
levelsikp21 = logspace(-2,2.1,50);
for kk = 1:length(ky)
    for nn = 1:length(n_new)
        if imag(kpara21(nn,kk))<=1.0e-2
            kpara21(nn,kk) = real(kpara21(nn,kk)) + 1.0e-2i;
        end
    end
end
contourf(log10(n_new),abs(ky),imag(kpara21)',levelsikp21,'Linecolor','none')
ic21=colorbar;
colormap(magma)
hold on
ic21.Ruler.TickLabelFormat = '%1.1f';
ic21.Ruler.Exponent = 1;
contour(log10(n_new),abs(ky),imag(kpara21)',[0.0013,0.0013],'-.w','Linewidth',3);
title('Im[k_{|| 21}] (m^{-1})')
set(gca,'ytick',[])
ylabel('k_{y} (m^{-1})')
xlabel('log_{10}n')
set(gca,'Fontsize',30,'FontName','CMU Serif')
set(gca,'Position',[0.55 0.12 0.3 0.35])
set(gca,'colorscale','log')
set(ic21,'Position',[h2(1) 0.12 h2(3) h2(4)])
text(0.02,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',30,...
    'Color','white')
hold off

%%

check = vpa(subs(wave_eq,kz,kz_arr(1,1)));

[V, D] = eig(check);

eig_V = V(:,2)./norm(V(:,2));

DIP = pi/2.0;

KEX = ky*eig_V(3,1)- kz_arr(1,1)*eig_V(2,1);
KEY = kz_arr(1,1)*eig_V(1,1)- kx*eig_V(3,1);
KEZ = kx*eig_V(2,1)-ky*eig_V(1,1);

% !	Calculate Field aligned component of Curl E
CurlE = KEX*cos(DIP) + KEZ*sin(DIP);
% !	Calculate Div(E)
DivE = kx*eig_V(1,1) + ky*eig_V(2,1) + kz_arr(1,1)*eig_V(3,1);

SX = (eig_V(2,1)*conj(KEZ)-eig_V(3,1)*conj(KEY))/(mu0*om);
SY = (eig_V(3,1)*conj(KEX)-eig_V(1,1)*conj(KEZ))/(mu0*om);
SZ = (eig_V(1,1)*conj(KEY)-eig_V(2,1)*conj(KEX))/(mu0*om);


%%
%--
% plot kx's

imme = find(imag(m)==0);
imme = imme(end)+1;


%%
%--
% transform data for log plot
y1 = sign(k1).*log10(abs(k1));
y2 = sign(k2).*log10(abs(k2));
y3 = sign(k3).*log10(abs(k3));
y4 = sign(k4).*log10(abs(k4));


%%

x0 = 0;
y0 = 0;
width = 1000;
height = 500;

figure(2)

set(gcf,'Position',[x0 y0 width height],'Color','w')
subplot(1,2,1)
semilogx(n_new, real(k1),'.k')
hold on
semilogx(n_new, real(k3),'.k')
semilogx(n_new, imag(k1),'.r')
semilogx(n_new, imag(k3),'.r')
ylim([-40 40])
xlabel('$n$ (m$^{-3}$)','Interpreter','latex')
ylabel('$k_{z}$ (m$^{-1}$)','Interpreter','latex')
set(gca,'Fontsize',25)

subplot(1,2,2)
semilogx(n_new, real(k2),'.k')
hold on
semilogx(n_new, real(k4),'.k')
semilogx(n_new, imag(k2),'.r')
semilogx(n_new, imag(k4),'.r')
ylim([-40 40])
xlabel('$n$ (m$^{-3}$)','Interpreter','latex')
%         ylabel('$k_{z}$ (m$^{-1}$)','Interpreter','latex')
yticks([])

set(gca,'Fontsize',25)

%%

figure(3)

set(gcf,'Position',[x0 y0 width height],'Color','w')
subplot(1,2,1)
plot(ax, real(k1),'.k')
hold on
plot(ax, real(k3),'.k')
plot(ax, imag(k1),'.r')
plot(ax, imag(k3),'.r')
ylim([-40 40])
xlabel('Position (m)','Interpreter','latex')
ylabel('$k_{z}$ (m$^{-1}$)','Interpreter','latex')
set(gca,'Fontsize',25)

subplot(1,2,2)
plot(ax, real(k2),'.k')
hold on
plot(ax, real(k4),'.k')
plot(ax, imag(k2),'.r')
plot(ax, imag(k4),'.r')
ylim([-40 40])
xlabel('Position (m)','Interpreter','latex')
%         ylabel('$k_{z}$ (m$^{-1}$)','Interpreter','latex')
yticks([])

set(gca,'Fontsize',25)

%%
figure(7)
plot(log10(n_new),real(y1),'.k')

hold on

plot(log10(n_new),imag(y1),'.r')
plot(log10(n_new),real(y3),'dk','MarkerSize',3)
plot(log10(n_new),imag(y3),'dr','MarkerSize',3)
legend('Re[k_{||}]', 'Im[k_{||}]', 'Re[k_{||}]', 'Im[k_{||}]')
xlabel('log_{10}|n |','Fontsize',16)
% vline(log10(n_new(imme)),'--k')
% yticklabels({'-10$^{3}$','-10$^{2}$','-10$^{1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'})
ylabel('log_{10}|k_{||}|','Fontsize',16)
set(gca,'XDir','reverse');

hold off

figure(8)
plot(log10(n_new),real(y2),'.k')

hold on

plot(log10(n_new),imag(y2),'.r')
plot(log10(n_new),real(y4),'dk','MarkerSize',3)
plot(log10(n_new),imag(y4),'dr','MarkerSize',3)
plot(log10(Nmax)*ones(1,npts),log10(linspace(1.0e-4,1.0e4,npts)),'b','Linewidth',3)
plot(log10(1.0e17)*ones(1,npts),log10(linspace(1.0e-4,1.0e4,npts)),'b','Linewidth',3)
legend('Re[k_{||1}]', 'Im[k_{||1}]', 'Re[k_{||1}]', 'Im[k_{||1}]')
xlabel('log_{10}|n |','Fontsize',16)
% vline(log10(n_new(imme)),'--k')
% yticklabels({'-10$^{3}$','-10$^{2}$','-10$^{1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'})
ylabel('log_{10}|k_{||}|','Fontsize',16)
set(gca,'XDir','reverse');

hold off

%%

figure(7)
plot(xax,real(y1),'.k')

hold on

plot(xax,imag(y1),'.r')
plot(xax,real(y3),'dk','MarkerSize',3)
plot(xax,imag(y3),'dr','MarkerSize',3)
legend('Re[k1]', 'Im[k1]', 'Re[k3]', 'Im[k3]')
xlabel('Position','Fontsize',16)
% vline(log10(N0(imme)),'--k')
% yticklabels({'-10$^{3}$','-10$^{2}$','-10$^{1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'})
ylabel('log_{10}|k_x|','Fontsize',16)

hold off

figure(8)
plot(xax,real(y2),'.k')

hold on

plot(xax,imag(y2),'.r')
plot(xax,real(y4),'dk','MarkerSize',3)
plot(xax,imag(y4),'dr','MarkerSize',3)
legend('Re[k2]', 'Im[k2]', 'Re[k4]', 'Im[k4]')
xlabel('Position','Fontsize',16)
% vline(log10(N0(imme)),'--k')
% yticklabels({'-10$^{3}$','-10$^{2}$','-10$^{1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'})
ylabel('log_{10}|k_x|','Fontsize',16)

hold off