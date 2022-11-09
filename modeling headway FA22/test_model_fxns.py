#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 15:55:59 2022

@author: vanessa
"""

# the purpose of this file is to test each piece of my model

#%%
"""
Unit tests for each piece of my model:
    
    calculation of a global timestep
    calculation of slope
    calculation of channel width
    calculation of depth
    calculation of bed exposure
    calculation of hydraulic radius transport efficiency
    calculation of sediment transport rate
    calculation of plucking rate
    calculation of abrasion rate
    calculation of attrition rate
    calculation of lateral sediment supply
    calculation of sedimentation rate
    calculation of total erosion rate
    
    """
    
# imports   
from numpy.testing import assert_allclose, assert_equal, assert_raises
import numpy as np


# test global timestep
def test_global_timestep():
    
    # parameters
    dx = 1000 # grid spacing, m
    x = np.arange(0, 2000, dx) # domain length, m
    x_node = x + (dx/2) # node-centered x, m
    kqs = 0.041 # sediment discharge coefficient
    r = 10. # runoff rate, m/yr
    area = (1./3.) * (x**2) # area set by Hack's law, m^2
    Q = r * area # total basin discharge, m^3/yr
    B = (2./3.) * x_node # basin width at any point along channel, m
    
    dt_global = 0.2 * (0.2 * dx * dx / (kqs*(Q[-1]/B[-1]))) # default timestep, yrs
    
    predicted_dt = 292.69

    assert_allclose(dt_global, predicted_dt, atol=1e-2)


# test slope
def test_slope():
    
    # parameters
    dx = 1000 # grid spacing, m
    x = np.arange(0, 2000, dx) # domain length, m
    eta = np.linspace(1, 0.1, len(x))
    
    slope = 0.0009
    
    predicted_slope = np.abs(np.diff(eta)/dx)

    assert_allclose(slope, predicted_slope, atol=1e-5)
    

# test channel width
def test_channel_width():
    
    # parameters
    dx = 1000 # grid spacing, m
    x = np.arange(0, 2000, dx) # domain length, m
    kb = 8.3e-8
    r = 10. # runoff rate, m/yr
    area = (1./3.) * (x**2) # area set by Hack's law, m^2
    Q = r * area # total basin discharge, m^3/yr
    S = 0.0009
    D = 0.03
    
    b = 0.015
    
    predicted_b = (kb * Q[1:] * (S ** (7/6))) / (D**(3/2))
    
    assert_allclose(b, predicted_b, atol=1e-3)
    
    
# test channel depth
def test_channel_depth():
    
    # parameters
    depth_coeff = 0.09801 # constant density and shear stress factor, eqn 9, W&S
    S = 0.0009
    D = 0.03
    
    depth = 3.267
    
    predicted_depth = (depth_coeff * D) / S
    
    assert_allclose(depth, predicted_depth, atol=1e-3)
    
    
# test bed exposure
def test_bed_exposure():
    
    # parameters
    Hstar = 0.5
    H = 200
    
    alpha = 0
    
    predicted_alpha = np.exp(-H/Hstar)
    
    assert_allclose(alpha, predicted_alpha, atol=1e-7)
    
    
# test hydraulic radius transport efficiency
def test_Rh_transport_efficiency():
    
    # parameters
    width = 0.015
    depth = 3.267
    
    Rh_efficiency = 0.00229
    
    predicted_Rh_efficiency = 1 - np.exp(-width/(2*depth))
    
    assert_allclose(Rh_efficiency, predicted_Rh_efficiency, atol=1e-5)
    
    
    
    
    
    
    
    
    
    