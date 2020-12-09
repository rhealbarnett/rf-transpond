# `rf-transpond` # 

A coupled 1D plasma transport (`transport_1d.m`) and cold plasma wave solver (`rf_wave_sol.m`) to determine ponderomotive effects in a fusion plasma close to radio-frequency antenna. 

## Dependencies ##

Matlab **[r2018b tested]**  

All subdirectories in this repository must be in the Matlab path.

## Running the verification tests ##

Both the transport and wave verification tests can be run via 
```
runtests('rf_transp_tests.m')
=======================================
Running wave solver dispersion verification
=======================================

...

=======================================
Running steady state transport MMS
=======================================

Iterations to tolerance information will be printed.

=======================================
Running time dependent transport MMS
=======================================

Max number of iterations, time step, total time taken information will be printed.
```

These tests will produce a figure for each verification test that will be saved in the `outputs/` subdirectory. 

### Wave solver verification ###

The input script for the wave solver verification is `wave_verification.m`. This will require the functions

- `dielec_tens.m`: calculates the cold plasma dielectric tensor,
- `plasma_freq.m` and `cyclo_freq.m`: calculate the plasma and cyclotron frequencies, respectively,
- `kz_dispersion.m`: calculates the cold plasma dispersion relation, 
- `fft_kz.m`: calculates the spatial FFT of the electric field solution from `rf_wave_sol.m`,
- `kz_spectrum.m`: calls `fft_kz.m` over a range of densities,

which are all housed in the `functions/` subdirectory. 
