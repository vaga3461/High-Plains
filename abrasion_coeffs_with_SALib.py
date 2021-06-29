#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 13:24:46 2021

@author: vanessa
"""


#%%
# this file will convert abrasion coefficient calibration efforts from 
# .ipynb to .py format
# it will also use SALib to explore parameter spaces

#%%

# import libraries
import numpy as np
import matplotlib.pyplot as plt
# from celluloid import Camera
import SALib
from SALib.sample import morris
from SALib.analyze import morris

#%%

# set up grid stuff
dx = 1
x = np.arange(0, 100, dx)
qs = np.zeros(len(x))
qs_ig = np.zeros(len(x))
qs_mtm = np.zeros(len(x))
qs_sed = np.zeros(len(x))

# define some constants
q = x
c = 1
porosity = 0.45
S = 0.001
H = 2

# lithologic percentages for a given basin
percent_ig = 0.58
percent_mtm = 0.41
percent_sed = 0.01

# set runtime
num_steps = 1000
dt = (0.2 * dx * dx)/c

#%%

# use SALib to generate conbinations of coefficients
problem = {
    'num vars': 3,
    'names': ['beta_ig', 'beta_mtm', 'beta_sed'],
    'bounds': [[2.5, 7.5],
               [7.6, 10],
               [11, 15]]}

abrasion_values = SALib.sample.morris.sample(problem, 5)

#%%

# write the model as a function because THE FUTURE IS MODULAR

def lith_predict(x,
                 q,
                 num_steps,
                 amt_ig,
                 amt_mtm,
                 amt_sed,
                 beta_ig, 
                 beta_mtm, 
                 beta_sed, 
                 c = 1, 
                 S = 0.001, 
                 porosity = 0.45, 
                 H = 2):
    
    """ The arguments to this function require the following data types:
        
        x : array
        q : array of length x
        num_steps : int
        amt_ig : float, < 1
        amt_mtm : float, < 1
        amt_sed : float, < 1
        (amt_ig + amt_mtm + amt_sed should always = 1)
        beta_ig : float
        beta_mtm : float
        beta_sed : float
        c : float or int
        S : float or int
        porosity : float, < 1
        H : float or int """
        
   # create arrays to hold sediment fluxes
    qs = np.zeros(len(x))
    qs_ig = np.zeros(len(x))
    qs_mtm = np.zeros(len(x))
    qs_sed = np.zeros(len(x))
    
    # create arrays to hold percentages of different lithologies
    theta = 1
    theta_ig = np.zeros(len(x))
    theta_mtm = np.zeros(len(x))
    theta_sed = np.zeros(len(x))
    theta_sed[:] = theta
    
    for i in range(num_steps):
        
        # calculate total sed flux and set boundary condition
        qs[1:] = c * q[1:] * S
        qs[0] = 0
    
        # calculate flux of each lithology
        qs_ig[1:] = qs[1:] * (theta_ig[:-1])
        qs_mtm[1:] = qs[1:] * (theta_mtm[:-1])
        qs_sed[1:] = qs[1:] * (theta_sed[:-1])
    
        # set constant, mixed feed of grains
        qs_ig[0] = amt_ig
        qs_mtm[0] = amt_mtm
        qs_sed[0] = amt_sed
    
        # update percentage of each grain type
        theta_ig[:-1] += ((-1/(porosity * H)) * ((np.diff(qs_ig)/dx) + (qs_ig[1:] * beta_ig))) * dt
        theta_mtm[:-1] += ((-1/(porosity * H)) * ((np.diff(qs_mtm)/dx) + (qs_mtm[1:] * beta_mtm))) * dt
        theta_sed[:-1] += ((-1/(porosity * H)) * ((np.diff(qs_sed)/dx) + (qs_sed[1:] * beta_sed))) * dt
    
        # conserve mass
        theta_total = theta_ig + theta_mtm + theta_sed
        # dtheta = theta - theta_total
        # theta_ig += dtheta * theta_ig
        # theta_mtm += dtheta * theta_mtm
        # theta_sed += dtheta * theta_sed
        
    # make a plot of the final output
    plt.plot(x, theta_ig/theta_total, label = 'igneous')
    plt.plot(x, theta_mtm/theta_total, label = 'metamorphic')
    plt.plot(x, theta_sed/theta_total, label = 'sedimentary')
    plt.legend()
    plt.show()
    
    return x, theta_ig, theta_mtm, theta_sed, theta_total
    

#%%

# now we loop over the generated values in our model

for i, X in enumerate(abrasion_values):
    
    # run the model with different parameter combos
    x, theta_ig, theta_mtm, theta_sed, theta_total = lith_predict(x, 
                                                                  q,
                                                                  num_steps,
                                                                  percent_ig, 
                                                                  percent_mtm, 
                                                                  percent_sed, 
                                                                  abrasion_values[i, 0], 
                                                                  abrasion_values[i, 1],
                                                                  abrasion_values[i, 2])
    
    
    
