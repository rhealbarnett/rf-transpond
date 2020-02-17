fslow = 10.0e6;
ffast = 10.0e6;

Afast = 1.0;
Aslow = 1.0;

Tfast = 1.0/ffast;
Tslow = 1.0/fslow;
tmax = 7*Tslow;
dt = Tfast/20;
tax = linspace(0,tmax,tmax/dt);

NP = floor(tmax/dt);

tseries_fast = zeros(1,NP);
tseries_slow = zeros(1,NP);

df = 1.0/(tmax);
fmax = 1.0/(2.0*dt);
freq_arr = linspace(1,fmax,fmax/df);


for ii=1:NP
    
    tseries_fast(1,ii) = Afast*cos(2.0*pi*ffast*dt*ii);
    tseries_slow(1,ii) = Aslow*sin(2.0*pi*fslow*dt*ii);
    
end

tseries_comb = tseries_fast.*tseries_slow;

fft_fast = fft(tseries_fast)/NP;
fft_slow = fft(tseries_slow)/NP;
fft_comb = fft(tseries_comb)/NP;

%%

% laplace transform test

xmax = 2.0*pi;
xmin = 0;
npts = 250;
t = linspace(xmin,xmax,npts);
wave_len = (xmax - xmin)/10.;
alpha = -0.5;
beta = 2.0*pi/wave_len;

% sig = exp(beta*t./100);
% sig = (1/2)*exp(1i*beta*t) + (1/2)*exp(-1i*beta*t);

for ii=1:npts
 
    signal(1,ii) = exp(alpha*t(1,ii) + 1i*beta*t(1,ii));
    
end

% f1 = @(t) exp((alpha + 1i*beta)*t);
% Lapf1 = @(s) LapTrans(f1,s);
% % t = linspace(xmin,xmax,npts);
% sc = 1:0.1:10;
% sp = 1:10;
% % yc = 2./sc.^3;
% yc = 1.0 ./ (sc - (alpha + 1i*beta));
% yp = arrayfun(Lapf1, sp);
% figure(1)
% plot(sc,yc,'r',sp,yp,'*');
% legend('Exact Transformation', 'Approximate Transformation');
% title('Numerical Laplace Transform (e^{(\alpha + i\beta) t})')
syms s

% signal(t) = exp(t);
% sig(t) = sin(beta*t);

f = exp(-s.*t);
test = trapz(t,(real(signal).*f));

% output = subs(test,s,beta);

% for ii=1:npts
%     
%     output(1,ii) = subs(test,s,signal(1,ii));
%     
% end

% lap_analytic = (s)./(s.^2 + beta^2);
lap_analytic = 1.0 / (s - beta);

s = linspace(2,10.0*pi,npts);
lap_analytic = subs(lap_analytic,s);
test = subs(test,s);

figure(1)
plot(s,test)
title('')
hold on

figure(1)
hold on
plot(s,lap_analytic)
legend('numerical laplace','analytical laplace')

%%

figure(2)
plot(real(output),imag(output),'.r','Markersize',10)
hold on
box on
xlim([-1.0 1.0])
ylim([-1.0 1.0])
xlabel('Real(s)')
ylabel('Imag(s)')
plot(linspace(-1,1,npts),zeros(1,npts),'--k')
plot(zeros(1,npts),linspace(-1,1,npts),'--k')



