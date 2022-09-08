#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:51:03 2022

@author: vanessa
"""

# this notebook contains the two-lithology fluvial model
# modified from the "basic" version to include
    # a formulation that honors basin geometry when calculating sediment flux
    
#%%
# LIBRARIES
import numpy as np
import matplotlib.pyplot as plt

#%%
# MODEL PARAMETERS

# gridding stuff (m)
dx = 1000
x = np.arange(0, 100000, dx)

# topographic information (m)
Hstar = 0.1 # characteristic sediment thickness
H = Hstar + np.zeros(len(x)) # sediment layer
etab = -H # bedrock elev
eta = etab + H # total elev

# abrasion coefficients from Attal and Lave 2006 in % per km
beta_ig = np.zeros(len(x))
beta_ig[:25] = 0.00004
beta_sed = np.zeros(len(x))
beta_sed[25:] = 0.00014
atr_factor = 0.00004

# erodibility values
k_ig = np.zeros(len(x))
k_ig[:25] = 0.0001
k_sed = np.zeros(len(x))
k_sed[25:] = 0.001

# width stuff
kB = 79.06  # From W&S
PB = 0.1  # from W&S
B = kB * x**PB

# discharge stuff
kh = 0.3
area = kh * x**2
runoff = 10.
Q = runoff * area

# runtime stuff
nsteps = 10000
dt = 1000

# uplift
U = 0.001

# transport and flow factors
I = 0.01
kQs = 0.041
#%%
Qs = np.zeros(len(x))
Q = np.zeros(len(x))
z = dx + np.zeros(len(x))
psi = 100.
porosity = 0.55
sedimentation_rate = np.zeros(len(x))

for i in range(nsteps):
    
    # apply uplift and set boundary conditions
    eta[-1] = U * dt
    etab[-1] = eta[-1]
    
    # calc slopes
    S = -np.diff(eta)/dx
    
    #  calculate e factor
    efac = np.exp(- H / Hstar)
    
    # calc Qs
    Qs[1:] = kQs * I * Q[:-1] * S**(7./6.) * (1.0 - efac[:-1])
    Qs[0] = 0
    dQsdx = np.diff(Qs) / dx
    
    # calc bedrock erosion
    #  calc bedrock erosion from stream power (plucking)
    ero_plucking_ig = efac[:-1] * (k_ig[1:] * Q[1:] * S)
    ero_plucking_sed = efac[:-1] * (k_sed[1:] * Q[1:] * S)
    
    #  calc bedrock erosion from abrasion
    ero_ab_ig = efac[:-1] * (beta_ig[:-1] * Qs[1:])   # <== change indexing: qs[1] represents node 0
    ero_ab_sed = efac[:-1] * (beta_sed[:-1] * Qs[1:])
    
    bedrock_ero = ero_plucking_ig + ero_plucking_sed + ero_ab_ig + ero_ab_sed
    
    #  calc change in bedrock elev
    etab[:-1] -= bedrock_ero * dt
    
    #  calc grain attrition rate
    atr = atr_factor * Qs[1:]
    
    # calc lateral sediment influx
    qL = z[:-1] - eta[:-1] / psi
    
    #  calc rate of change in alluvial thickness
    sedimentation_rate[:-1] = -((1 / porosity) * (dQsdx + atr - ero_plucking_ig - qL))
    
    #  update sediment thickness
    H[:-1] += sedimentation_rate[:-1] * dt
    H[H < 0] = 0
    
    #  update elev
    eta[:-1] = etab[:-1] + H[:-1]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    