#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 15:54:05 2022

@author: vanessa
"""

# this code attempts to write a model 
# only slightly more complex than the one Greg already wrote
# it will basically just take his code and try to add in a single bedrock erosion process
# assume all one lithology

#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
# Code to try this
kh = 1./3.  # Hack coefficient
psi = 100.0  # response-time coefficient, years per meter
U = 0.0001  # baselevel lowering rate, m/y
nsteps = 10000
dt = 1000.0  # time step, years
kB = 79.06  # From W&S
PB = 0.1  # from W&S
runoff = 10.0  # bankfull runoff coefficient, m/yr
dx = 1000.0  # node spacing, m
I = 0.01  # intermittency factor
kQs = 0.041  # transport coefficient

x = dx * np.linspace(1, 10, 10)  # steamwise distance, m

A = kh * x**2  # drainage area, m2
Q = runoff * A  # discharge, cmy
B = kB * x**PB  # valley width, m

z = dx + np.zeros(len(x))  # stream height, m
H = dx + np.zeros(len(x))  # adjacent topography avg. height, m
Qs = np.zeros(len(x))

char_sed = 0.1
sed = char_sed + np.zeros(len(x))

k = 0.001

bedrock = -sed

elev = bedrock + sed

for i in range(nsteps):
    
    # lower baselevel
    z[-1] -= U * dt
    
    # calc slopes
    S = -np.diff(z) / dx
    
    # calc bed exposure
    efac = np.exp(-sed / char_sed)
    
    # calc sed transport
    Qs[1:] = kQs * I * Q[:-1] * S**(7/6) * (1.0 - efac[:-1])
    
    # calc sed flux
    dQsdx = np.diff(Qs) / dx
    
    # calc bedrock erosion
    bedrock_ero = k * Q[1:] * S * efac[:-1]
    
    # calc lateral influx
    qL = (H[:-1] - z[:-1]) / psi
    
    # calc change in bed elev (sed thickness) from lateral influx
    dHdt = (qL - dQsdx) / B[:-1]
    
    
    elev[:-1] += delev_dt * dt # update bed elevation
    dHdt = -qL / (2 * kh * x[:-1])  # rate of sed loss from adjacent topo, m/yr
    H[:-1] += dHdt * dt  # update avg elevation of adjacent topo
    
plt.plot(x, z, '.-')
plt.plot(x[:-1], H[:-1], '.-')
plt.xlabel('Distance (m)')
plt.ylabel('z, H (m)')
plt.legend(['Stream height', 'Adjacent basin height'])