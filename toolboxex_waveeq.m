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

m = 1.0;
c = 1.0;
a = 0.0;
f = 0.0;

numberofPDE = 1;
model = createpde(numberofPDE);

sq1 = [3 4 -.8 .8 .8 -.8  .8 .8 -.8 -.8]';
gdm = [sq1];
g = decsg(sq1);
geometryFromEdges(model,g);

figure
pdegplot(model,'EdgeLabels','on');
axis([-.9 .9 -.9 .9]);
title 'Block Geometry With Edge Labels Displayed'



