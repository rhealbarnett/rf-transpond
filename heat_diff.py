#-------------------------------------------#
# solve heat diffusion eqn                  #
# eqn 5.82 in Callen text, chap 5 pg 27     #
# explicit and implicit methods             #
# rlb 170117                                #
#-------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import math 
import itertools

colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
cc = itertools.cycle(colours)

#%%

#--- spatial domain
# domain min and max in metres
xmin = 0.
xmax = 1.0   

# number of points 
nptsx = 21   

# spatial axis
xax = np.linspace(xmin, xmax, nptsx) 
 
# spatial step size (dx) of 0.05 m
dx = 1.0 / (nptsx-1)  

#--- temporal bits
#- courant (CFS) condition -#
# generally, I'm used to a courant condition defined in terms of EM wave propagation
# rather than heat (ie c*dt/dx <= 1.0/sqrt(n) for c the vacuum light speed and n the 
# number of spatial dimensions). I realise this isn't exactly applicable for heat 
# diffusion. With the help of professor google and Numerical Recipes (the art of
# scientific computing) chapter 17, I've managed to find something more useful;
# (2*D)*(dt/(dx)^2) <= 1.0 where D is (2/3)*chi for eq 5.82 from the Callen text

# diffusion coefficient-- arbitrarily chose 1 for the moment (for simplicity)
chi = 1.

D = (2.0 / 3.0)*chi

# time step -- explicit dt needs to be set using the Courant condition
dt = 0.999*((dx**2) / (2.0*D))

# simulation time (seconds)-- kept small, otherwise the number of iterations is quite large
tmin = 0.
tmax = 1.0

# number of time steps
nptst = int(tmax / dt)

# time axis
tax = np.linspace(tmin, tmax, nptst)

#--- calculate coefficients, initialise arrays, boundary & initial conds
# coefficient on rhs
alpha = D*(dt / (dx**2))

# temperature array
temp = np.zeros(nptsx)

# boundary conds - chose what seems to be a commonly used boundary conditions 
# to allow testing
T_left = 0.
T_right = 0.

# initial conds - as above for the boundary conds
temp = np.sin(np.pi*xax/xmax)

temp[0] = T_left
temp[nptsx-1] = T_right

## exact solution for error analysis
#exact = np.sin(np.pi*xax/xmax)*np.exp(-1.0*D*np.pi**2*tax/(xmax**2))

# initial min and max temps
temp_max = max(temp)
temp_min = min(temp)

plt.figure(1)
plt.plot(xax[0:nptsx], temp, 'k--', label='t=0.00s')
plt.hold(True)

#--- solver 
for ii in range(1,nptst):
    
    temp[1:nptsx-1] = temp[1:nptsx-1] + alpha*(temp[2:nptsx] - 2.0*temp[1:nptsx-1] + temp[0:nptsx-2])
    
    #plot every 100 iterations
    if math.fmod(int(ii), 100)==0:
        c = next(cc)
	plt.plot(xax[0:nptsx], temp, color=c, label='t=%0.2fs' % (ii*dt))
        plt.hold(True)
        
plt.legend(loc='upper right',prop={'size':10})    
plt.xlabel('x')
plt.ylabel('Temp')
plt.ylim(temp_min, temp_max)
plt.xlim(xmin, xmax)    
plt.title('Explicit (FTCS) Method')
plt.show(1)

#%%

#------------------------------------------------------------#

# --- implicit (BTCS) method --- #

# redefine number of temporal steps -- no longer needs to obey the CFS condition so can be arbitrary
nptst = 21

# new dt (again not defined by CFS)
dt = 1.0/(nptst-1)

# new time axes      
tax = np.linspace(tmin, tmax, nptst)

# calculate alpha again based on new dt value       
alpha = D*(dt / (dx**2))

# initialise coefficient matrix for Ax = b
A = np.zeros((nptsx, nptsx))

# fill tridiagonal matrix 
for ii in range(1,nptsx-1):
#    A[ii,ii-1] = -1.0*r
#    A[ii,ii] = 1.0 + 2.0*r
#    A[ii,ii+1] = -1.0*r
    A[ii,ii-1] = -1.0*alpha
    A[ii,ii] = 1.0 + 2.0*alpha
    A[ii,ii+1] = -1.0*alpha

# boundary conditions T1 = T_left and Tn = T_right
A[0,0] = 1.0
A[-1,-1] = 1.0

# reset initial conditions
temp = np.sin(np.pi*xax/xmax)
temp[0] = T_left
temp[nptsx-1] = T_right

# plot initial function
plt.figure(2)
plt.plot(xax[0:nptsx], temp, 'k--', label='t=0.00s')
plt.hold(True)

#--- solver
for ii in range(1,nptst):
    
    t = ii*dt
    print 'time %0.2f s' % (ii*dt)
    
    # set up b column vector as the 'old' temperature values
    b = temp
    
    # calculate new temperature values
    temp_new = np.linalg.solve(A, b)
    
    # fill old temperature array with new values for next iteration
    temp = temp_new
    
    # plot every 0.25 seconds
    if math.fmod(int(ii), 5)==0:
        c = next(cc)
        plt.plot(xax[0:nptsx], temp, color=c, label='t=%0.2fs' % (ii*dt))
        plt.hold(True)
        
plt.legend(loc='upper right',prop={'size':10})    
plt.xlabel('x')
plt.ylabel('Temp')
plt.ylim(temp_min, temp_max)
plt.xlim(xmin, xmax) 
plt.title('Implicit (BTCS) Method') 
plt.show(2)  
    
#%%
    
#------------------------------------------------------------#
# --- implicit (Crank-Nicholson) method --- #

# initialise coefficient matrix and column vector for Ax = b
A = np.zeros((nptsx, nptsx))
b = np.zeros(nptsx)

# set initial conditions again
temp = np.sin(np.pi*xax/xmax)
temp[0] = T_left
temp[nptsx-1] = T_right

# plot initial conds
plt.figure(3)
plt.plot(xax[0:nptsx], temp, 'k--', label='t=0.00s')
plt.hold(True)

for ii in range(1,nptsx-1):
#    A[ii,ii-1] = -1.0*r
#    A[ii,ii] = 1.0 + 2.0*r
#    A[ii,ii+1] = -1.0*r
    A[ii,ii-1] = -1.0*alpha
    A[ii,ii] = 2.0*(1.0 + alpha)
    A[ii,ii+1] = -1.0*alpha
# boundary conditions T1 = T_left and Tn = T_right
A[0,0] = 1.0
A[-1,-1] = 1.0

#--- solver
for ii in range(1,nptst):
    
    t = ii*dt
    print 'time %0.2f s' % (ii*dt)
    
    # calculate b column vector values
    b[1:nptsx-1] = alpha*temp[2:nptsx] + 2.0*(1.0-alpha)*temp[1:nptsx-1] + alpha*temp[0:nptsx-2]
    
    # calculate new temperature values
    temp_new = np.linalg.solve(A, b)
    
    # fill old temperature array with new values for next iteration
    temp = temp_new
    
    # plot every 0.25 seconds
    if math.fmod(int(ii), 5)==0:
        c = next(cc)
        plt.plot(xax[0:nptsx], temp, color=c, label='t=%0.2fs' % (ii*dt))
        plt.hold(True)
        
plt.legend(loc='upper right',prop={'size':10})    
plt.xlabel('x')
plt.ylabel('Temp')
plt.ylim(temp_min, temp_max)
plt.xlim(xmin, xmax)    
plt.title('Implicit (Crank-Nicholson) Method')
plt.show(3)
    

    
    
    
    
    
    
    
    
    
 
