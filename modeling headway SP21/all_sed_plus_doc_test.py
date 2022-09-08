#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 15:36:28 2022

@author: vanessa
"""

#%%
# this file contains a special version of my model
# for the case where sediment supply is infinite (H >> H*)
# and attrition does not occur (beta = 0)

#%%
# start by importing libraries
import numpy as np
import matplotlib.pyplot as plt

#%%
# now set up arrays and parameters
dx = 1000 # grid spacing
x = np.arange(0, 3000, dx) # domain length

# H = 100 + np.zeros(len(x)) # sediment thickness
# z = np.linspace(1, 0.1, len(x)) # + H # need to start with slight bedrock slope
z = np.zeros(len(x))

U = 0.0005 # uplift rate
phi = 0.55 # sediment porosity
kqs = 0.041 # sediment discharge coefficient
I = 0.01 # intermittency factor
r = 10. # runoff rate
kxb = 25 # valley width coeffecient
Pxb = 1/5 # valley width exponent
# Hstar = 0.1 # characteristic sediment thickness

B = kxb * (x**Pxb) # valley width 
Q = B * r * x # total discharge

#%%
# now create a function that attempts to capture this all as simply as possible
# don't need bedrock slope because that's only for calculating channel width, 
# and we're not even worrying about that as a comparison metric

def all_sed(dx, x, z, U, phi, kqs, I, r, kxb, Pxb, B, Q, num_steps=550000):
    
    """This is a simplified version of my fluvial model. Below is a test.

    Examples
    --------
    >>> test = [dx, x, z, U, phi, kqs, I, r, kxb, Pxb, B, Q]
    >>> S, Qs, E, dzdt, z, dt, total_dzdt = all_sed(*test)
    >>> S[1]
    0.0831
    """
    
    # set timestep
    dt = (0.5 * dx * dx / (kqs*Q[-1]))
    
    # create arrays
    Qs = np.zeros(len(x))
    E = np.zeros(len(x))
    dzdt = np.zeros(len(x))
    
    # set boundary conditions
    Qs[0] = 0
    E[-1] = 0
    dzdt[-1] = 0
    
    # track uplift
    total_dzdt = 0
    
    for i in range(num_steps):
        
        # calculate slope
        S = np.abs(np.diff(z)/dx)
        
        # calculate sediment transport
        Qs[1:] = kqs * I * Q[1:] * S**(7./6.)
        
        # erosion
        E[:-1] = ((1/((1-phi)*B[1:])) * np.diff(Qs)/dx)
        
        # calculate change in elevation
        dzdt[:-1] = U - E[:-1]
        
        # update profile
        z[:] += dzdt * dt
        
        # track total change in elev so we can account for this when comparing to prediction
        total_dzdt += dzdt * dt
        
    cum_time = num_steps * dt
    # print(cum_time)
    # print(S[1])
    
        
    return (S, Qs, E, dzdt, z, dt, total_dzdt)

#%%
# now try writing a doc test
"""This is a simplified version of my fluvial model. Below is a test.

Examples
--------
>>> dx = 1000
>>> x = np.arange(0, 10000, dx)
>>> z = np.zeros(len(x))
>>> U = 0.0005
>>> phi = 0.55
>>> kqs = 0.041
>>> I = 0.01
>>> r = 10.
>>> kxb = 25
>>> Pxb = 1/5
>>> B = kxb * (x**Pxb) # valley width 
>>> Q = B * r * x # total discharge
>>> test = [dx, x, z, U, phi, kqs, I, r, kxb, Pxb, B, Q]
>>> S, Qs, E, dzdt, z, dt, total_dzdt = all_sed(*test)
>>> S[:]
0.0831
"""
