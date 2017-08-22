# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 19:31:52 2016

@author: luisalejo
"""

import numpy as np
from matplotlib import pyplot

T = 100.
dt = 0.02
N = int(T/dt) + 1
t = np.linspace(0.0, T, N)

z0 = 100.
b0 = 10.
zt = 100.
g = 9.81
u = np.array([z0, b0])
z = np.zeros(N)
z[0] = z0
b = np.zeros(N)
b[0] = b0

for n in range(1, N):
    u = u + dt*np.array([u[1], g*(1-u[0]/zt)])
    z[n] = u[0]
    b[n] = u[1]

pyplot.figure(figsize=(10,4))   #set plot size
pyplot.ylim(40,160)             #y-axis plot limits
pyplot.tick_params(axis='both', labelsize=14) #increase font size for ticks
pyplot.xlabel('t (s)', fontsize=14) #x label
pyplot.ylabel('z (m)', fontsize=14) #y label
pyplot.plot(t,z, 'r-');

pyplot.figure(figsize=(10,4))   #set plot size
pyplot.ylim(-15,15)             #y-axis plot limits
pyplot.tick_params(axis='both', labelsize=14) #increase font size for ticks
pyplot.xlabel('t (s)', fontsize=14) #x label
pyplot.ylabel('b (m/s)', fontsize=14) #y label
pyplot.plot(t,b, 'g-');


# Definining the different values for dt
dt_values = np.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0001])

# Array containing the z solutions for each dt (for each grid).
z_values = np.empty_like(dt_values, dtype=np.ndarray)

for i, dt in enumerate(dt_values):
    N = int(T/dt)+1  #Number of time-steps
    t = np.linspace(0.0, T, N)
    
    # initila conditions
    u = np.array([z0, b0])
    z = np.empty_like(t)
    z[0] = z0
    
    # time loop
    for n in range(1,N):
        u = u +dt*np.array([u[1], g*(1-u[0]/zt)])
        z[n] = u[0]
        
    z_values[i] = z.copy()
