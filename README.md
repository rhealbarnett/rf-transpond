# `rf-transpond` # 

A coupled 1D plasma transport (`transport_1d.m`) and cold plasma wave solver (`rf_wave_sol.m`) to determine ponderomotive effects in a fusion plasma close to radio-frequency antenna. 

## Dependencies ##

Matlab **[r2018b tested]**
All subdirectories in this repository must be in the Matlab path.

## Running the verification tests ##

Both the transport and wave verification tests can be run via 

`runtests('rf_transp_tests.m')`
`=======================================`
`Running wave solver dispersion verification`
`=======================================`
`...`
`=======================================`
`Running steady state transport MMS`
`=======================================`
`
`**Iterations to tolerance information will be printed**`
`
`=======================================`
`Running time dependent transport MMS`
`=======================================`
`
`**Max number of iterations, time step, total time taken information will be printed**`

