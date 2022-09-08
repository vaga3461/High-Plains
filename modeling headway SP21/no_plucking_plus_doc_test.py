#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 13:16:57 2022

@author: vanessa
"""

#%%
# this file contains a special version of my model
# for the case where bedrock erosion occurs only via abrasion, no plucking
# this means new coarse sediment cannot be generated in the model

#%%
# start by importing libraries
import numpy as np
import matplotlib.pyplot as plt

#%%
# define parameters

dx = 1000 # grid spacing
x = np.arange(0, 10000, dx) # domain length
r = 10. # runoff rate
kh = 1/3 # hack coefficient
h = 2 # hack exponent
Hstar = 0.5 # characteristic sediment thickness
# beta = np.zeros(len(x))
# beta[:3] = 0.0004
# beta[3:] = 0.004
beta = 0.004 # abrasion coefficient
# psi = 0.00004 # attrition factor
# K = np.zeros(len(x))
# K[:3] = 0.000001
# K[3:] = 0.00001
# K = 0.0001 # rock erodibility
gamma = 0.5 # fraction of coarse sediment from hillslope erosion
kxb = 25 # valley width coefficient
Pxb = 1/5 # valley width exponent
kb = 8.3e-8 # channel width coefficient
phi = 0.55 # sediment porosity
D = 0.05 # grain size
U = 0.0005

#%%
# set up arrays

H = Hstar + np.zeros(len(x)) # sediment thickness
etab = np.zeros(len(x)) # bedrock elevation array
etab[:] = np.linspace(1, 0.1, len(x)) # need to start with slight bedrock slope
eta = etab + H # total topographic elevation

#%%
# calculate constant, derivative values

B = kxb * (x**Pxb) # valley width
# Q = r * kh * (x**h) # total discharge
Q = r * B * x

# define more constants
kqs = 0.041 # sediment discharge coefficient
I = 0.01 # intermittency factor

#%%
numerator = (U*gamma*x)*(1 - np.exp(-(beta * x)/3))
denominator = beta*(1 + (beta*x))*(kqs*I*Q)

slope_pred = (numerator/denominator)**(6./7.)

#%%
def no_plucking(x,
          dx,
          Hstar,
          H,
          etab,
          eta,
          beta,
          gamma,
          kxb,
          Pxb,
          kb,
          phi,
          D,
          U,
          Q,
          B,
          num_steps = 7000000):
    
    # define more constants
    kqs = 0.041 # sediment discharge coefficient
    I = 0.01 # intermittency factor
    
    # calculate timestep
    dt_global = (0.2 * dx * dx / (kqs*Q[-1]))
    run_duration = dt_global * num_steps  # <== here's how long we want to run
    cum_time = 0.0  # <== keep track of elapsed time
    
    # define arrays
    b = np.zeros(len(x)) # channel width
    Eb_a = np.zeros(len(x)) # abrasion rate
    Eb = np.zeros(len(x)) # bedrock erosion rate
    Eh = np.zeros(len(x)) # sedimentation rate
    E = np.zeros(len(x)) # total erosion rate
    q = np.zeros(len(x)) # unit discharge
    Qs = np.zeros(len(x)) # total sediment transport
    qs = np.zeros(len(x)) # unit sediment transport
    ql = np.zeros(len(x)) # lateral sediment supply
    
    # set boundary conditions
    b[0] = 0
    H[-1] = 0
    E[-1] = U
    q[0] = 0
    Qs[0] = 0
    qs[0] = 0
    ql[0] = 0
    
    while cum_time < run_duration:  # <== use a while loop because dt varies by iteration
        
        # first calculate rates
            
        # calculate slope
        S = np.abs(np.diff(eta)/dx)
        
        # calculate channel width (L)
        b = (kb * Q[1:] * (S ** (7 / 6))) / (D**(3/2))
        
        # calculate unit discharge (L^2/T)
        q[1:] = Q[1:]/b
        
        # calculate bed exposure
        alpha = np.exp(-H/Hstar)
        
        # calculate sediment transport (L^3/T and L^2/T)
        # Qs[1:] = c * I * Q[1:] * S**(7/6) * (1 - alpha[:-1])
        Qs[1:] = kqs * I * Q[1:] * np.sign(S) * (np.abs(S)) ** (7/6) * (1-alpha[:-1])
        qs[1:] = Qs[1:]/b
        
        # calculate individual erosion mechanism rates (L^2/T)
        Eb_a[:-1] = beta * Qs[1:] * alpha[:-1]
        
        # calculate total bedrock erosion rate (L/T)
        # Eb[:-1] = Eb_p[:-1] + Eb_a[:-1]
        Eb[:-1] = Eb_a[:-1] / B[1:]
        
        # calculate attrition rate (L^2/T)
        atr = beta * Qs
        
        # calculate lateral sediment inputs (L^2/T)
        ql = ((Eb * gamma) / beta) * (1 - np.exp(-beta * (x/3)))
        
        # calculate sedimentation rate (L/T)
        Eh[:-1] = - (1/ ((1 - phi) * B[1:])) * ((np.diff(Qs)/dx) + atr[1:] - ql[:-1])
        #Eh[1:] = - (1/ (1 - phi) * B[1:]) * ((np.diff(Qs)/dx) + atr[1:] - Eb_p[:-1] - ql[:-1])
        
        # calculate total erosion rate (L/T)
        E[:-1] = Eb[:-1] + Eh[:-1]
        
        
        
        # Calculate maximum allowable time-step size
        
        #  set adaptive timestep
        #  first check time to flat surface
        elev_diff = np.diff(eta)/dx
        ero_diff = np.diff(E)/dx
        #valid_places = np.where(ero_diff < 0)
        valid_places = np.where(ero_diff < 0)[0]  # <== we just want the array, not the full tuple from where()
        if len(valid_places) > 0:  # <== in case there ARE no locations...
            times_to_flat = np.abs(elev_diff[valid_places]/ero_diff[valid_places])
        else:
            times_to_flat = np.array([dt_global])  # <== ...we just revert to the global dt
        min_time_to_flat = np.amin(times_to_flat)

        #  then check time to deplete all sediment
        #sed_depletion_locations = np.where(sedimentation_rate < 0)
        sed_depletion_locations = np.where(Eh < 0)[0]  # <== we just want the array, not the full tuple from where()
        if len(sed_depletion_locations) > 0:  # <== in case there ARE no locations...
            times_to_no_sed = np.abs(H[sed_depletion_locations]/Eh[sed_depletion_locations])
        else:
            times_to_no_sed = np.array([dt_global])  # <== ...we just revert to the global dt
        min_time_to_no_sed = np.amin(times_to_no_sed)

        #  check for smaller condition
        dt = min(min_time_to_flat, min_time_to_no_sed)

        #  if larger than global step size, limit to global
        dt = min(dt, dt_global)
        
        
        
        # Update quantities
        
        # update boundary conditions
        eta[-1] -= U * dt
        etab[-1] = eta[-1]
        
        # update topography
        etab[:-1] -= Eb[:-1] * dt
        H[:-1] += Eh[:-1] * dt
        H[H<0] = 0
        eta[:-1] = etab[:-1] + H[:-1]
        
        # Advance time
        cum_time += dt
        
        # print(dt, "timestep in years")
        
        if any(E[:] != U):
            continue
        else:
            break
            
    print(cum_time, "years")
        
    return (S, b, q, Qs, qs, Eb, atr, ql, Eh, E, eta, etab, H, dt)

#%%
# now try writing a doc test
"""This is a simplified version of my fluvial model. Below is a test.

Examples
--------
>>> dx = 1000
>>> x = np.arange(0, 10000, dx)
>>> H = Hstar + np.zeros(len(x))
>>> etab = np.zeros(len(x))
>>> etab[:] = np.linspace(1, 0.1, len(x)) 
>>> eta = etab + H 
>>> U = 0.001
>>> gamma = 0.5
>>> beta = 0.004
>>> kqs = 0.041
>>> I = 0.01
>>> r = 10.
>>> kxb = 25
>>> Pxb = 1./5.
>>> kb = 8.3e-8 
>>> phi = 0.55 
>>> D = 0.05
>>> B_4 = kxb * (x[4]**Pxb)
>>> Q_4 = r * B_4 * x[4]
>>> test = [x, dx, Hstar, H, etab, eta, beta, gamma, kxb, Pxb, kb, phi, D, U, Q_4, B_4]
>>> S, b, q, Qs, qs, Eb, atr, ql, Eh, E, eta, etab, H, dt = no_plucking(*test)
>>> S[4]
0.0251
"""