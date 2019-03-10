%------------------------------------------------------------%
% example of particle position/trajectory in electric field
% simple ponderomotive 
%------------------------------------------------------------%

const = constants();
me = const.me;
e = const.e;

me = me/me;
e = e/e;

tmin = 0;
tmax = 1.0;
t = linspace(tmin,tmax,1000);
npts = length(t);
dt = (tmax - tmin)/(npts-1);

period = tmax/10.0;
freq = 1.0/period;
om = 2.0*pi*freq;
Eamp = ones(1,npts);
E = Eamp.*cos(om*t);

a = e*E/me;
v = dt*a;
x = dt*v;

for ii=1:npts/5
    figure(1)
    plot(x(ii),0,'.b','Markersize',20);
    legend(['time = ' num2str(double(ii)*dt/period) ' T'])
    xlim([min(x) max(x)])
    set(gca,'YTickLabel',[])
    pause(0.05)
end

xax = linspace(0,1,npts);
dx = (max(xax) - min(xax))/(npts-1);
Eamp = linspace(0,0.5,npts);
E = cos(om*t) + Eamp;
% E = ((E.^2));
% Ediff = (E(2:npts) - E(1:npts-1))/(dx);

% a =  - (e^2/(4.0*me^2*om^2))*Ediff;
a = e*E/me;
v = dt*a;
x = dt*v;

% x = zeros(1,npts);

% x(1) = 0;
% 
% x(2) = x(1) + v(1)*dt;
% 
% for ii=2:npts-1
%     x(ii+1) = dt^2*(e*E(ii+1)/me) + 2.0*x(ii) - x(ii-1);
% %     v(ii) = (x(ii+1) - x(ii-1))/(2.0*dt);
% end

for ii=1:npts
    
    figure(1)
    plot(x(ii),0,'b.','Markersize',20);
    legend(['time = ' num2str(double(ii)*dt/period) ' T'])
    xlim([-max(x) max(x)])
    set(gca,'YTickLabel',[])
    ylim([-1,1])
    txt = {'Low E Amplitude \rightarrow'};
    text(0,0.5,txt,'Fontsize',14)
    pause(0.05)

end


