% ------------------------------------------------------------------- %
% Run example of coupling RF wave and transport solvers.
% Use LAPD as example device. 
% 201126 rlbarnett
% ------------------------------------------------------------------- %

function ans = lapd_example()
    
    lapd_equil
    wave_sol_plots
    plots = 0;
    transport_1d
    wave_sol_plots
    
end

