#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:17:23 2022

@author: vanessa
"""

#%%
"""
Unit tests for special cases 1 and 2:
    
    All (unlimited) sediment
    Case 1: No abrasion or attrition; Case 2: sediment attrition
    No bedrock
    
    """
    
# imports   
from numpy.testing import assert_allclose, assert_equal, assert_raises
import numpy as np

# define the model for SC1
def test_sc1():
    
    # parameters
    dx = 1000 # grid spacing
    U = 0.001 # uplift rate
    phi = 0.55 # sediment porosity
    kqs = 0.041 # sediment discharge coefficient
    I = 0.01 # intermittency factor
    r = 10. # runoff rate
    kxb = 25 # valley width coeffecient (a, above)
    Pxb = (1/5) # valley width exponent
    
    x = np.arange(0, 3000, dx) # domain length
    x_node = x + dx/2
    z = np.linspace(1, 0, len(x)) # need to start with slight bedrock slope
    
    B = kxb * (x_node**Pxb) # valley width   ## TWEAKED
    Q = (r * kxb * x**(6/5))/(1 + Pxb) # discharge  ## WE TWEAKED THIS
    
    # set timestep
    dt = (0.5 * dx * dx / (kqs*Q[-1]))
    
    # create arrays
    Qs = np.zeros(len(x))
    E = np.zeros(len(x))
    dzdt = np.zeros(len(x))
    
    # set boundary conditions
    Qs[0] = 0
    dzdt[-1] = 0 # dzdt at outlet = 0 for case where dzdt = 0
    
    num_steps = 6000000
    
    for i in range(num_steps):
        
        # calculate slope
        actual_slope_SC1 = np.abs(np.diff(z)/dx)
        
        # calculate sediment transport
        ##Qs[1:] = kqs * I * Q[1:] * S**(7./6.)
        Qs[1:] = kqs * I * Q[1:] * actual_slope_SC1**(7./6.)
        
        # erosion
        # "old" indexing from [:-1]
        E[:-1] = (1/((1-phi)*B[:-1])) * ((np.diff(Qs)/dx))
        
        # calculate rate of elevation change
        # "old" indexing of E from [:-1]
        dzdt[:-1] = U - E[:-1]
        
        # update profile
        z += dzdt * dt
        
    return (actual_slope_SC1, Qs, E, dzdt, z, dt)
    
    # analytical solution for slope
    predicted_slope_SC1= ((U * (1-phi) / (kqs * I * r))**(6./7.))
    
    assert_allclose(actual_slope_SC1[1], predicted_slope_SC1, rtol = 1e-5, atol=0)
    
    print("predicted slopes SC1: ", predicted_slope_SC1)
    print("actual slopes SC1: ", actual_slope_SC1)
    
    
    
    
    
# define the model for SC2
def test_sc2():

    # parameters
    dx = 1000 # grid spacing
    U = 0.001 # uplift rate
    phi = 0.55 # sediment porosity
    kqs = 0.041 # sediment discharge coefficient
    I = 0.01 # intermittency factor
    r = 10. # runoff rate
    beta = 0.4 # abrasion coeffcient
    kxb = 25 # valley width coeffecient (a, above)
    Pxb = (1/5) # valley width exponent
    
    x = np.arange(0, 3000, dx) # domain length
    x_node = x + dx/2
    z = np.linspace(1, 0, len(x)) # need to start with slight bedrock slope
    
    B = kxb * (x_node**Pxb) # valley width   ## TWEAKED
    Q = (r * kxb * x**(6/5))/(1 + Pxb) # discharge  ## WE TWEAKED THIS
    
    # set timestep
    dt = (0.5 * dx * dx / (kqs*Q[-1]))
    
    # create arrays
    Qs = np.zeros(len(x))
    E = np.zeros(len(x))
    dzdt = np.zeros(len(x))
    
    # set boundary conditions
    Qs[0] = 0
    dzdt[-1] = 0 # dzdt at outlet = 0 for case where dzdt = 0
    
    num_steps = 6000000
    
    for i in range(num_steps):
        
        # calculate slope
        actual_slope_SC2 = np.abs(np.diff(z)/dx)
        
        # calculate sediment transport
        ##Qs[1:] = kqs * I * Q[1:] * S**(7./6.)
        Qs[1:] = kqs * I * Q[1:] * actual_slope_SC2**(7./6.)
        
        # erosion
        # "old" indexing from [:-1]
        E[:-1] = (1/((1-phi)*B[:-1])) * ((np.diff(Qs)/dx) + (beta * Qs[1:]))
        
        # calculate rate of elevation change
        # "old" indexing of E from [:-1]
        dzdt[:-1] = U - E[:-1]
        
        # update profile
        z += dzdt * dt
        
    return (actual_slope_SC2, Qs, E, dzdt, z, dt)
    
    # analytical solution for slope
    predicted_slope_SC2 = ((U * B[1:] * x[1:] * (1 - phi))/(1 + (beta * x[1:])) * 1/(kqs*I*Q[1:]))**(6/7)
    
    assert_allclose(actual_slope_SC2[1], predicted_slope_SC2, rtol = 1e-5, atol=0)
    
    print("predicted slopes SC2: ", predicted_slope_SC2)
    print("actual slopes SC2: ", actual_slope_SC2)



