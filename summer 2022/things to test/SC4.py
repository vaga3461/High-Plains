#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 15:31:57 2022

@author: vanessa
"""

#%%
# start by importing libraries
import numpy as np
import matplotlib.pyplot as plt

#%%
# define parameters

dx = 1000 # grid spacing
x = np.arange(0, 5000, dx) # domain length
x_node = x + dx/2

baselevel_rate = 0.001 # uplift rate
phi = 0.55 # sediment porosity
kqs = 0.041 # sediment discharge coefficient
I = 0.01 # intermittency factor
r = 10. # runoff rate
S0 = 0.15 # slope linearizer
K = 0.00001 # erodibility
# K[:3] = 0.000001
# K[3:] = 0.00001
Hstar = 0.5
kxb = 25 # valley width coeffecient (a, above)
Pxb = (1/5) # valley width exponent

#%%
# set up arrays

H = Hstar + np.zeros(len(x)) # sediment thickness
etab = np.zeros(len(x)) # bedrock elevation array
etab[:] = np.linspace(1, 0.1, len(x)) # need to start with slight bedrock slope
eta = etab + H # total topographic elevation

#%%
# calculate constant, derivative values

B = kxb * (x_node**Pxb) # valley width   ## TWEAKED
Q = (r * kxb * x**(6/5))/(1 + Pxb) # discharge  ## WE TWEAKED THIS

#%%
def model(x,
          dx,
          Hstar,
          H,
          etab,
          eta,
          kqs,
          I,
          K,
          kxb,
          Pxb,
          phi,
          baselevel_rate,
          Q,
          B,
          num_steps = 200000):
    
#     dx = 1000 # grid spacing
#     x = np.arange(0, 5000, dx) # domain length
#     x_node = x + dx/2

    dx = 1000 # grid spacing
    x = np.arange(0, 5000, dx) # domain length
    x_node = x + dx/2
        
    # calculate timestep
    dt = 0.2 * (0.2 * dx * dx / (kqs*(Q[-1]/B[-1])))
    
    # define arrays
    Eb = np.zeros(len(x)) # bedrock erosion rate
    Eb_total = np.zeros(len(x)) # total bedrock erosion rate
    Eh = np.zeros(len(x)) # sedimentation rate
    E = np.zeros(len(x)) # total erosion rate
    Qs = np.zeros(len(x)) # total sediment transport
    
    # set boundary conditions
    H[-1] = 0
    E[-1] = 0 # also try baselevel_rate
    Qs[0] = 0
    
    for i in range(num_steps):
            
        # calculate slope
        S = np.abs(np.diff(eta)/dx)
        
        # calculate bed exposure
        alpha = np.exp(-H/Hstar)
        
        # calculate sediment transport (L^3/T)
        Qs[1:] = kqs * I * Q[1:] * (S** (7./6.)) * (1-alpha[:-1])
        
        # calculate bedrock erosion rate (plucking) (L^2/T)
        Eb[:-1] = K * Q[1:] * S * alpha[:-1]
        
        # total bedrock erosion (L/T)
        Eb_total = Eb / B
       
        # calculate sedimentation rate (L/T)
        Eh[:-1] = - (1/((1 - phi) * B[1:])) * ((np.diff(Qs)/dx) - Eb[:-1])
        
        # calculate total erosion rate (L/T)
        E[:-1] = Eb_total[:-1] + Eh[:-1]
        
        # update boundary conditions
        eta[-1] -= baselevel_rate * dt
        etab[-1] = eta[-1]
        
        # update topography
        etab[:-1] -= Eb[:-1] * dt
        H[:-1] += Eh[:-1] * dt
        H[H<0] = 0
        eta[:-1] = etab[:-1] + H[:-1]

        
    return (S, Qs, Eb, Eb_total, Eh, E, eta, etab, H, dt)

#%%
# now write a test
test = [dx, x, Hstar, H, etab, eta, kqs, I, K, kxb, Pxb, phi, baselevel_rate, Q, B]

# and run it
S, Qs, Eb, Eb_total, Eh, E, eta, etab, H, dt = model(*test)

#%%
E

#%%
plt.plot(x[1:], S)