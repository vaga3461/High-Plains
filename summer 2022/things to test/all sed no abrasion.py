#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 10:36:42 2022

@author: vanessa
"""

# the purpose of this exercise is to get the constant slope test for SC1 working
# this means changing B to be calculated at nodes rather than faces
# and setting up a unit test. 
# (Alternatively, might set up a unit test that allows us to set a tolerance)
# we'll also use this space to make some figures for a visual representation of the accuracy of our analytical solution

# as a reminder, SC1 us all (unlimited) sediment, no abrasion/attrition, no bedrock erosion

#%%
# start by important libraries
import numpy as np
import matplotlib.pyplot as plt

#%%
# set up arrays and parameters

dx = 1000 # grid spacing
x = np.arange(0, 3000, dx) # domain length
x_node = x + dx/2

z = np.linspace(1, 0, len(x)) # need to start with slight bedrock slope

U = 0.001 # uplift rate
phi = 0.55 # sediment porosity
kqs = 0.041 # sediment discharge coefficient
I = 0.01 # intermittency factor
r = 10. # runoff rate
kxb = 25 # valley width coeffecient (a, above)
Pxb = (1/5) # valley width exponent

B = kxb * (x_node**Pxb) # valley width   ## TWEAKED
Q = (r * kxb * x**(6/5))/(1 + Pxb) # discharge  ## WE TWEAKED THIS

#%%
# now create a function for the case where dzdt=0

def all_sed(dx, x, z, U, phi, kqs, I, r, kxb, Pxb, B, Q, num_steps=6000000):
    
    """This is a simplified version of my fluvial model. Below is a test.

    Examples
    --------
    >>> test = [dx, x, z, U, phi, kqs, I, r, kxb, Pxb, B, Q]
    >>> S, Qs, E, dzdt, model_z, dt = all_sed(*test)
    >>> S[:]
    0.15
    """
    
    # set timestep
    dt = (0.5 * dx * dx / (kqs*Q[-1]))
    
    # create arrays
    Qs = np.zeros(len(x))
    E = np.zeros(len(x))
    dzdt = np.zeros(len(x))
    
    # set boundary conditions
    Qs[0] = 0
    dzdt[-1] = 0 # dzdt at outlet = 0 for case where dzdt = 0
    
    for i in range(num_steps):
        
        # calculate slope
        S = np.abs(np.diff(z)/dx)
        
        # calculate sediment transport
        ##Qs[1:] = kqs * I * Q[1:] * S**(7./6.)
        Qs[1:] = kqs * I * Q[1:] * S**(7./6.)
        
        # erosion
        # "old" indexing from [:-1]
        E[:-1] = (1/((1-phi)*B[:-1])) * ((np.diff(Qs)/dx))
        
        # calculate rate of elevation change
        # "old" indexing of E from [:-1]
        dzdt[:-1] = U - E[:-1]
        
        # update profile
        z += dzdt * dt
        
    cum_time = num_steps * dt
    # print(cum_time)
        
    return (S, Qs, E, dzdt, z, dt)

#%%
# now write a test
test = [dx, x, z, U, phi, kqs, I, r, kxb, Pxb, B, Q]

# and run it
S, Qs, E, dzdt, model_z, dt = all_sed(*test)

#%%
# predicted slope calculation and line

slope = ((U * (1-phi) / (kqs * I * r))**(6./7.))
line = (-slope * x) + model_z[0]

#%%
# make a comparison plot

plt.plot(x, model_z)
plt.plot(x, line)

























