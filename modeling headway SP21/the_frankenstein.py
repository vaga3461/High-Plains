#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:29:22 2021

@author: vanessa
"""


#%%
# this .py file will frankenstein together the stream power and abrasion model
# with the lithology prediction model

#%%
# start by importing libraries
import numpy as np
import matplotlib.pyplot as plt
import SALib
from SALib.sample import morris
from SALib.analyze import morris

#%%
