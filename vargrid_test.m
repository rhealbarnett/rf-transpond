%---------------------------------------------------------------------%
% A program to test 1D grid refinement.
% adapted from page 2 of 
% https://www3.nd.edu/~gtryggva/CFD-Course/2011-Lecture-25.pdf
% rlbarnett c3149416 160119
%---------------------------------------------------------------------%

% sets xmax
L = 0.1; 

% for internal layer calculation
% sign change at xc
xc = L/2.0; 

% "strength": concentrates points either at the centre of the domain (A>0)
% or at the boundaries of the domain (A<0)
A = 2.5;

% root order
ro = 0.5;

% initialise xi parameter array
% spacing in the centre currently is 0.5*dx
smax = 1.0;
smin = 0.0;
% s = smin:2*2.4e-4:smax; 
s = linspace(smin,smax,2048-1);

% length of the x parameter, i.e. the number of points
n = length(s);

% calculate the x values from xi
x = L*(s.^(1/ro));
% x = L*s+A*(xc-L*s).*s.*(1-s);

% plots x as a funtion of xi
figure(1);
plot(s,x,'LineWidth',2);
hold on

% loop plots the lines along x(xi), visually shows dx sizing effectively
for ii=1:round(n/50):n
    
    plot([s(ii),s(ii)],[0.0, x(ii)]);
    plot([0.0,s(ii)],[x(ii), x(ii)]);
    
end

% plot parameters/labels
xlabel('\xi','Fontsize',18) 
ylabel('x','Fontsize',18)
set(gca,'Box','on'); 
set(gca,'Fontsize',18, 'LineWidth',2)
plot([smin,smax],[0,L],'r')
% text(0.05,L-0.1,['A=',num2str(A),...
% ' x=',num2str(xc)],'Fontsize',18)

hold off;
