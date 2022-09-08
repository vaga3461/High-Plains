#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:26:43 2022

@author: vanessa
"""

# the purpose of this exercise is to get the constant slope test for SC2 working
# this means changing B to be calculated at nodes rather than faces
# and setting up a unit test. 
# (Alternatively, might set up a unit test that allows us to set a tolerance)
# we'll also use this space to make some figures for a visual representation of the accuracy of our analytical solution

# as a reminder, SC2 us all (unlimited) sediment, with abrasion/attrition of grains, but no bedrock erosion

#%%
# start by important libraries
import numpy as np
import matplotlib.pyplot as plt

#%%
# now set up arrays and parameters
dx = 1000 # grid spacing

U = 0.001 # uplift rate
phi = 0.55 # sediment porosity
kqs = 0.041 # sediment discharge coefficient
I = 0.01 # intermittency factor
r = 10. # runoff rate
beta = 0.0004 # abrasion factor
kxb = 25 # valley width coeffecient
Pxb = 1/5 # valley width exponent


x = np.arange(0, 10000, dx) # domain length
x_node = x + dx/2
z = np.linspace(1, 0, len(x)) # need to start with slight bedrock slope
 
B = kxb * (x_node**Pxb) # valley width   ## TWEAKED
Q = (r * kxb * x**(6/5))/(1 + Pxb) # discharge  ## WE TWEAKED THIS

#%%
# now define a function

def all_sed_attrition(dx, x, z, U, phi, kqs, I, r, kxb, Pxb, beta, B, Q, num_steps=7000000):
    
    # set timestep
    dt = (0.5 * dx * dx / (kqs*Q[-1]))
    
    # create arrays
    Qs = np.zeros(len(x))
    E = np.zeros(len(x))
    dzdt = np.zeros(len(x))
    
    # set boundary conditions
    Qs[0] = 0
    #E[-1] = 0
    dzdt[-1] = 0

    for i in range(num_steps):
    
        # calculate slope
        S = np.abs(np.diff(z)/dx)
    
        # calculate sediment transport
        ##Qs[1:] = kqs * I * Q[1:] * S**(7./6.)
        Qs[1:] = kqs * I * Q[1:] * S**(7./6.)
    
        # erosion
        # "old" indexing from [:-1]
        E[:-1] = (1/((1-phi)*B[:-1])) * ((np.diff(Qs)/dx) + (beta * Qs[1:]))
    
        # calculate rate of elevation change
        # "old" indexing of E from [:-1]
        dzdt[:-1] = U - E[:-1]
    
        # update profile
        z += dzdt * dt
    
    return (S, Qs, E, dzdt, z, dt)

#%%
# and design a test

test = [dx, x, z, U, phi, kqs, I, r, kxb, Pxb, beta, B, Q]

# and run it
slope, Qs, E, dzdt, model_z, dt = all_sed_attrition(*test)

#%%
# calc the slope with the analytical soln

slope_pred = ((U * B[1:] * x[1:] * (1 - phi))/(1 + (beta * x[1:])) * 1/(kqs*I*Q[1:]))**(6/7)

#%%
print(slope)
print(slope_pred)

print(slope - slope_pred)

#%%
plt.plot(x[1:], slope)
plt.plot(x[1:], slope_pred)
plt.show()

#%%
topo_pred = -(slope_pred * x[1:]) + model_z[0]

#%%
plt.plot(x[1:], model_z[1:])
plt.plot(x[1:], topo_pred)
plt.show()


















