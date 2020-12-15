import numpy as np
import matplotlib.pyplot as plt
# this script only knows the down_to_east function if you import it
from scripts.functions import down_to_east

class erosion_model(object):
    
    def __init__(self,
                 channel_length = 1000,
                 erodibility = 0.00001,
                 slope = 0.002,
                 zmax = 1000,
                 m = 0.5,
                 n = 1,
                 dx = 20,
                 uplift = 0.001,
                 duration = 1E6,
                 timestep = 25):
        
        self.channel_length = channel_length
        self.erodibility = erodibility
        self.x = np.arange(0, channel_length, dx)
        self.uplift = down_to_east(self.x, 0.003)
        self.area = 0.5*self.x
        self.zmax = zmax
        self.slope = slope
        self.z = (slope*channel_length) + zmax
        self.m = m
        self.n = n
        self.duration = duration
        self.timestep = timestep
        
        print('it works!')
        
    def run_one_step(self):
        # this doesn't seem to work because z is a single value so diff fails
        dzdx = np.diff(self.z)/dx
        absolute = np.abs(dzdx)
        self.erosion = self.erodibility*(self.area[1:]**self.m)*(dzdx**self.n)
        dzdt = self.uplift - self.erosion
        self.z[1:] += dzdt * self.timestep
        self.z[0] += self.zmax + (self.uplift * self.timestep)
        self.z[self.z<0] = 0
        
    def run_n_steps(self, n):
        for i in range(n):
            self.run_one_step()
    
    def run_model(self):
        print("testing the 'run_model' function")
        print("member variable z: %.1f" % self.z)
            
#     def run_model(self,
#                   channel_length = 1000,
#                   erodibility = 0.00001,
#                   slope = 0.002,
#                   zmax = 1000,
#                   m = 0.5,
#                   n = 1,
#                   dx = 20,
#                   uplift = 0.001,
#                   duration = 1000000,
#                   timestep = 25):
    
#         model = erosion_model(channel_length = channel_length,
#                               erodibility = erodibility,
#                               slope = slope,
#                               zmax = zmax,
#                               m = m,
#                               n = n,
#                               dx = dx,
#                               uplift = uplift,
#                               duration = duration,
#                               timestep = timestep)
    
#         def init():
#             ax.set_ylim(0, 2 * np.amax(model.z))
#             ax.set_xlim(0, model.channel_length)
#             ax.set_ylabel('Height (m)')
#             ax.set_xlabel('Distance (m)')
#             return(obj)
    
#         fig, ax = plt.subplots()
#         xdata = []
#         ydata = []

#         def update(i):
#             ax.cla()
#             model.run_n_steps(save_every)
#             xdata = model.x
#             ydata = model.z
#             ax.plot(xdata, ydata)

#         def report_values(model):
#             """Report values of various quantities."""
#             print('total erosion: ' + str(model.erosion)) 

#         def make_plot(model):
#             plt.plot(self.x, self.z)
#             plt.show