# ** `rf-transpond` ** # 

#### vels.py

###### RK4 and Boris push method solvers to reproduce fig 13 in VE paper

Both methods: x and v results are an order of magnitude lower than expected??  
Shape of vx, vy and vz look correct though.  
Note: both the rk4 and Boris push method are second order accurate.

#### rk4 test.py

###### python's in-built rk4 method

Sort of. Not super easy to use.  
Works fine for a simple case, extension into 3-component position and velocity not working yet.  
Don't really understand what each step is doing just yet.  

#### disp rels.py

###### solves the dispersion relation

Fine for specific cases (eg right-handed wave = R). General n^2 equation with specified angle doesn't though? Can't figure out why. Double checked expressions for S, D, P, R, and L as well as A, B and C in Swanson textbook and they all agree.  
But the first term in A doesn't zero completely (nan in the array), which adds to the second term. Might be part of the problem. 

#### heat diff.py

###### solves 1d heat diffusion

Explicit, as well as implicit (BTCS and crank-nicholson) methods.  
Still to do:  check error with changing dt, solution around sharp features, 2d.  
