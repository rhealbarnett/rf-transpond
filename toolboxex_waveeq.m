%-------------------------------------%
% pde toolbox example                 %
% wave eqn on a square domain         %
% rlb c3149416 230218                 %
%-------------------------------------%

%----------------------------------------------------------------------%
%
%   wave equation form 
%   partial^2 u/partial t^2 - del^2 u = 0
%    
%   need to express in a form that the toolbox can interpret
%   ie of the form
%   m*partial^2 u/partial t^2 - del.(c*del u) + a*u = f
%   which means m = 1, c = 1, a = 0, f = 0
%
%----------------------------------------------------------------------%

%--
% specify number of pdes: don't think there is a limit to this
% create pde model
numberofPDE = 1;
model = createpde(numberofPDE);

%--
% define geometry: this can be done via command line/script, but can also
% be done visually in the pde Modeler toolbox app, and the parameters
% exported as gd, sf and ns (defined in the Matlab documentation 
% https://www.mathworks.com/help/pde/ug/decsg.html)
sq1 = [3 4 -.8 .8 .8 -.8  .8 .8 -.8 -.8]';
gdm = [sq1];
g = decsg(sq1);
geometryFromEdges(model,g);

%--
% create plot of geometry
figure(1)
pdegplot(model,'EdgeLabels','on');
axis([-.9 .9 -.9 .9]);
title 'Block Geometry With Edge Labels Displayed'

%--
% specify coefficients for the pde
% not 100% sure what the 'd' coefficient is
m = 1.0;
c = 1.0;
a = 0.0;
f = 0.0;

specifyCoefficients(model,'m',m,'d',0,'c',c,'a',a,'f',f);

%--
% specify the boundary conditions: can be either neumann or dirichlet
% conditions, where h*u=r is the dirichlet condition and n*c*grad(u) + q*u
% = g is the neumann condition
uLeftRight = applyBoundaryCondition(model,'dirichlet','Edge',[2,4],'u',0);
uTopBottom = applyBoundaryCondition(model,'neumann','Edge',([1 3]),'g',0);

%--
% generate the mesh: can set the min/max size of the finite elements
% can "jiggle" the mesh, which adjusts the centre of the element, but
% apparently is no longer recommended
hmax = 0.1;
generateMesh(model,'Hmax',hmax);

%--
% create plot of geometry with mesh
figure(2)
pdemesh(model);
axis([-.9 .9 -.9 .9]);
axis equal
xlabel x
ylabel y

%--
% set initial conditions: not quite sure why there are two initial
% condition functions... maybe to set different initial conditions at
% different locations?
u0 = @(location) atan(cos(pi/2*location.x));
ut0 = @(location) 3*sin(pi*location.x).*exp(sin(pi/2*location.y));
setInitialConditions(model,u0,ut0);

%--
% generate time array
n = 31;
tlist = linspace(0,5,n);

%--
% statistics will report:
% # successful steps
% # failed attempts
% # function evaluations
% # partial derivatives
% # LU decompositions
% # solutions of linear systems
model.SolverOptions.ReportStatistics ='on';

%--
% solve pde: not sure why they are taking the nodal solution, or what other
% solutions are possible
% Ah, printing the result shows that there is a nodal solution, solution
% times, (x,y,z)gradients and the mesh.
result = solvepde(model,tlist);
u = result.NodalSolution;

%--
% generates figure for each time step, plays like a movie -- pretty cool
figure(3)
umax = max(max(u));
umin = min(min(u));
for i = 1:n
    pdeplot(model,'XYData',u(:,i),'ZData',u(:,i),'ZStyle','continuous',...
                  'Mesh','off','XYGrid','on','ColorBar','off');
    axis([-1 1 -1 1 umin umax]);
    caxis([umin umax]);
    xlabel x
    ylabel y
    zlabel u
    M(i) = getframe;
end
