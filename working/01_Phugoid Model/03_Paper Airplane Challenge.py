# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 18:58:59 2017

@author: Alejo
"""

from math import sin, cos, log, ceil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import time
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16
start_time = time.time() 

#Model Parameters:
g = 9.8         #gravity
v_t = 4.9        #trim velocity
LtoD = 5.       #Lift to drag ratio
DtoL = 1./LtoD   #Drag to lift ratio
T = 100.
dt = 0.01
N = int(T/dt) + 1
samples = 20

#Set initial conditions
v0_values = np.linspace(3., 6.5, samples)
theta0_values = np.linspace(-np.pi/2, np.pi/2, samples)
x0 = 0.
y0 = 2.


#Definition of the function f(u)
def f(u):
    """Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : array of floats
        array containing the solution at time t.
        
    Returns
    -------
    dudt : array of floats
        array containing the RHS of f(u), given u.
    """
    
    v = u[0]
    theta = u[1]
    x = u[2]
    y = u[3]
    return np.array([-g*sin(theta) - DtoL*g/v_t**2*v**2,
                      -g*cos(theta)/v + g/v_t**2*v,
                      v*cos(theta),
                      v*sin(theta)])


#Definition of the Euler-step function
def euler_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of floats
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equations.
    dt : float
        time-increment.
    
    Returns
    -------
    u_{n+1} : array of floats
              approximate solution at the next time step.
    """
    
    return u + dt * f(u)

u_set = []
u_set_trimmed = []

#%%

#Start solving
for i, v0 in enumerate(v0_values):
    for j, theta0 in enumerate(theta0_values):
        
        #initialise the array containing the instantaneous solutions
        u = np.zeros((N,4))
        u[0] = np.array([v0, theta0, x0, y0])
        
        #time loop
        for n in range(N-1):
            u[n+1] = euler_step(u[n], f, dt)
            
            #Limiting the loop to values y>=0
            if u[n+1][3] < 0:
                u[n+1] = np.zeros((1,4))
                break
        #store the value of u for the given v0 and theta0
        u_set.append(u)


#Get rid of the zeros in the u array
for e in range(len(v0_values)*len(theta0_values)):
    mask = np.all(np.equal(u_set[e], 0), axis = 1)
    u_set_trimmed.append(u_set[e][~mask])


#Find out the largets x value    
x_max = 0.      #Initialise a max value
index_maxj = 0  #Initialise the index i where x_max is inside u_set_trimmed
index_maxi = 0  #Initialise the index j where x_max is inside U_set_trimmed[i]

#Iterate through the the values of x
for i, item in enumerate(u_set_trimmed):
    for j in range(len(item)):
        x_test = item[j][2]
        #Logic test to find x_max and its indexes
        if x_test > x_max:
            x_max = x_test
            index_maxj = j
            index_maxi = i


#Gather results
v0_opt = u_set_trimmed[index_maxi][0][0]
theta0_opt = u_set_trimmed[index_maxi][0][1]

#Show results and time
print('The maximum distance is ' + str('%.2f' % x_max) + ', achieved with an initial velocity and pitch angle of ' + str('%.2f' % v0_opt) + ' m/s and ' + str('%.2f' % theta0_opt) + ' rads.')

print("--- %s seconds ---" % (time.time() - start_time))

plt.figure(1)
plt.figure(figsize=(10, 4))
plt.ylim(0,4)
plt.xlim(-1,20)
plt.xlabel('x (m)', fontsize = 14)
plt.ylabel('y (m)', fontsize = 14)
plt.plot(u_set_trimmed[index_maxi][:,2], u_set_trimmed[index_maxi][:,3], 'k-')