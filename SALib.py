# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 11:44:44 2020

@author: vanes
"""

#%%
# starting to play with sensitivity analysis, using terrainbento and SALib
# going to start by testing only a few values of three parameters:
# low k, high k, v_s
# and evolving a very simple, small landscape with terrainbento BasicHy using SALib-generated
# parameter combinations. Then will analyze parameter space both by inspection of outputs
# (for now, just total material removed),
# and with SALib analysis tools (mean and variance of effect of parameters)

#%%
# import libraries

import SALib
from SALib.sample import morris
from SALib.analyze import morris

import numpy as np
np.random.seed(42)

import matplotlib
import matplotlib.pyplot as plt

import holoviews as hv

from terrainbento import BasicHy, Clock, NotCoreNodeBaselevelHandler

from landlab import imshow_grid, RasterModelGrid
from landlab.io.netcdf import write_netcdf
from landlab.values import random

#%%
# create ranges for the inputs

problem = {
    'num_vars': 3,
    'names': ['low_k', 'high_k', 'v_s'],
    'bounds': [[1.0e-7, 1.0e-5],
               [1.0e-6, 1.0e-4],
               [1.0e-2, 1.0e2]]
}

#%%
# generate samples
# using Morris sampler, where resulting matrix will have (G + 1) * T rows and D columns
# D is number of parameters (3)
# G is number of groups (if no groups selected, then just number of parameters again) (3)
# T is number of trajectories - start with 3
# so output should have 12 rows and 3 columns

param_values = SALib.sample.morris.sample(problem, 3)

# col 0 is low k
# col 1 is high k
# col 2 is v_s

#%%
# next step will be iterating through these combos (rows of matrix) within model
# notation like this:

# Y = np.zeros([param_values.shape[0]])

# for i, X in enumerate(param_values):
#     Y[i] = evaluate_model(X)

# where Y is one output that you want to measure (for example: total material removed)
# it has the same length of the matrix because you get one output for each row of the matrix
# in loop, X is one combo of parameters (so, one row) in param_values, and
# i iterates through those rows
# so Y[i] - evaluate_model(X) means the output of run i (1, 2, 3...) when evaluated with 
# parameters of row i

#%%
# looping over values from SALib

for i, X in enumerate(param_values):
    
    # create grid and elevation field
    grid = RasterModelGrid((25, 40), xy_spacing=40) 
    z = grid.add_zeros('node', 'topographic__elevation', units = 'm')
    z += -grid.node_x * 0.0045 + 4000 + np.random.rand(len(grid.node_y))
    z_original = z.copy()
    
    # set boundaries
    grid.set_closed_boundaries_at_grid_edges(bottom_is_closed = True,
                                             left_is_closed = True,
                                             right_is_closed = False,
                                             top_is_closed = True)

    # runtime and timestep of each model run
    clock = Clock(start=0, step=1000, stop=1e7)
    
    # set erodibility values
    erodibility = np.zeros(len(z))
    erodibility_grid = np.reshape(erodibility, (25, 40))
    erodibility_grid[:, 0:15] = param_values[i, 0] # low erodibility
    erodibility_grid[:, 15:] = param_values[i, 1] # high erodibility

    # instantiate the model with i combination of X 
    # (given combo of param_values at row i of param_values matrix)
    model = BasicHy(clock,
                    grid,
                    water_erodibility = erodibility_grid,
                    settling_velocity = param_values[i, 2])
    
    # run the model
    model.run()
    
    # track elevation change and store in numpy array, 
    # where each column is elevation for each model run
    elevation_change_tracker = np.empty((1000, 12))
    np.append(elevation_change_tracker[:, i], z_original - z)
    
    
    # look at some outputs
    plt.figure()
    imshow_grid(grid, 'topographic__elevation')
    
