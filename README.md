# `rf-transpond` # 

A coupled 1D plasma transport (`transport_1d.m`) and cold plasma wave solver (`rf_wave_sol.m`) to determine ponderomotive effects in a fusion plasma close to radio-frequency antenna. 

## Quickstart

```
git clone https://github.com/rhealbarnett/rf-transpond.git
cd rf-transpond
addpath(genpath('./'))
runtests
```

## `rf_wave_sol.m` ##

Solves the frequency domain cold plasma wave equation along one axis. The electric field solution is calculated along *z* using a second order central finite difference scheme, with wavenumbers required for *x* and *y*. The three equations corresponding to each electric field component are written in matrix form, and a solution found using Matlab's \ operator. 

## `transport_1d.m` ##

Solves the single fluid continuity and momentum equations in 1D (domain parallel to a confining magnetic field B0) on a staggered, non-uniform grid. Advection terms are discretised using a first order upwind scheme, second order terms are discretised with a second order central difference scheme. The time stepping is implicit-explicit (IMEX).    

## Dependencies ##

Matlab **[r2018b tested]**  

All subdirectories in this repository must be in the Matlab path.

## Running the verification tests ##

Both the transport and wave verification tests can be run via 
```
runtests('rf_transp_test.m')
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

Required parameters for the wave verification test are defined in the input script `wave_verification.m`. Running the verification via `runtests('rf_transp_tests.m')` requires the functions

- `dielec_tens.m`: calculates the cold plasma dielectric tensor,
- `plasma_freq.m` and `cyclo_freq.m`: calculate the plasma and cyclotron frequencies, respectively,
- `kz_dispersion.m`: calculates the cold plasma dispersion relation, 
- `fft_kz.m`: calculates the spatial FFT of the electric field solution from `rf_wave_sol.m`,
- `kz_spectrum.m`: calls `fft_kz.m` over a range of densities,

which are all housed in the `functions/` subdirectory, which must be in the Matlab path in order for the functions to be accessed. 

### Transport solver verification ###

Required parameters for both the steady state and time dependent method of manufactured (MMS) verification tests are defined in the input script `transport_mms.m`. Running the MMS verifications reqiures the functions

- `mms_source_cont.m`: calculates the MMS source term for the continuity equation,
- `mms_source_mom.m`: calculates the MMS source term for the momentum equation,
- `run_mms.m`: calls the transport script `transport_1d.m` for each set step size (spatial and temporal) and stores the calculated L2, L inifinity errors. 

which are also all housed in the `functions/` subdirectory.     

# Running the coupled wave / transport case for LAPD-like parameters #

Parameters for the transport and wave codes are defined in the input scripts `lapd_equil.m` and `lapd_rfparams.m` respectively. The equilibrium density and velocity profiles used as initial conditions are provided in `/inputs/equil_transport_input.mat`, which is loaded in `lapd_equil.m`. 

The test case can be run via 

```
lapd_example()
```

which calls `lapd_equil.m` and the transport solver script `transport_1d.m`, as well as the plotting routine, which saves an initial and final figure to `outputs/`. Running this example requires the functions

- `pressure_source_stag.m`: calculates the grad pressure term on the right hand side of the velocity update equation,
- `pond_source.m`: calculates the ponderomotive term on the right hand side of the velocity equation. In this specific example, the ponderomotive force is calculated using parallel gradients in all components of the electric field solution from `rf_wave_sol.m`, but this can be changed (details are contained in the function file).

as well as the previously listed `dielec_tens.m`, `plasma_freq.m` and `cyclo_freq.m`. The wave solver function `rf_wave_sol.m` is called inside the time stepping loop of `transport_1d.m`, providing the coupling between the plasma transport and cold plasma wave solvers.  
 







