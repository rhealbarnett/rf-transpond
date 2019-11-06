% ----------------------------------------------------------------------- %
% Script to read and plot matlab (.mat) files for the coupled 1D 
% wave / transport solver. 
%
% rlbarnett c3149416 191106
% ----------------------------------------------------------------------- %

filepath = '/Volumes/DATA/LAPD/matlab/results_jsource_kyzero/';

lapd_equib;

for ii=212:106:47382
    
    filename = strcat(filepath, 'coupled_transport_', num2str(ii),'.mat');
    
    transport.freq = freq;
    transport.period = period;
    transport.B0 = B0;
    transport.source_dist = source_dist;
    transport.scale_fact = scale_fact;
    transport.ey_source = ey_source;
    transport.R = R;
    transport.kx = kx;
    transport.ky = ky;
    
    save(filename,'-struct','transport','-append');
    
end
