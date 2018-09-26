
constants;
% transport_test;

transport_1d('staggered','central',...
  'neumann',0,'neumann',0,...
  'dirichlet',cs,'dirichlet',cs/2,const,'transport_test.m');