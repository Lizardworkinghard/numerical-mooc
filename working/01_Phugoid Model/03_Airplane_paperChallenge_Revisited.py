# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 11:05:41 2018

@author: villegas
"""


"""
Paper Airplane Challenge - Revisited
An attempt to apporach the problem from a diffecet perspective, trying to understand better the way the data
is stored in lists and ndarrays, as well as the best practices to slice these data sets and write scripst
in function of their lengths for a more general/flexible numerical model.
"""

import time
start_time = time.clock()
import numpy as np
import matplotlib.pyplot as plt



#Constants and discretization of the time domain
g = 9.81
v_t = 4.9
LtoD = 5.
DtoL = 1/LtoD
t = 100
dt = 0.01
N = int(t/dt) + 1

#Initial values. Ranges of values for v0 and theta0
v0_values = np.linspace(3.,6.5, num=2, endpoint=True)
theta0_values = np.linspace(-np.pi/2, np.pi/2, num=500, endpoint=True)
x0 = 0
y0 = 2.


#Definition of the function f as u'=f(u)
def f(u):
    v = u[0]
    theta = u[1]
    x = u[2]
    y = u[3]
    return np.array([-g*np.sin(theta)-DtoL*g*v**2/(v_t**2),
                     -g*np.cos(theta)/v + g*v/(v_t**2),
                     v*np.cos(theta),
                     v*np.sin(theta)])

    
#Defining the Euler method as u_n + dt*u'
def euler_step(u, f, dt):
    return u + dt*f(u)







#ITERATION-PROCESSING ----------------------------------------------------------------------------------
    
#initialising and empty LIST to hold all solutions
u_set = []

#Numerical process to get u_set

for i, v0 in enumerate(v0_values):
    for j, theta0 in enumerate(theta0_values):
        
        u = np.zeros((N,4))
        u[0] = np.array([v0, theta0, x0, y0])
       
        #Iterating through the range of initial values and breaking if y<0
        for n in range(N-1):
            u[n+1] = euler_step(u[n], f, dt)
            
            if u[n+1,3] < 0:
                u[n+1] = np.zeros((1,4))
                break
    
        u_set.append(u)



#ANALYSIS-POSTPROCESSING----------------------------------------------------------------------------------

#Finding the maximum value reached for x (u_set[i][j,2]) and recording its indexes
x_max = 0
x_max_index = (0,0)

for i in range(len(u_set)):
    for j in range(N-1):
        x = u_set[i][j,2]
        if x >= x_max:
            x_max = x
            x_max_index = (i, j)
        

#Extracting all values of x and y off the general solution u_set
x_values = np.zeros((len(u_set),N))
y_values = np.zeros((len(u_set),N))

for i in range(len(u_set)):
    for j in range(N-1):
        x = u_set[i][:,2]
        y = u_set[i][:,3]
    
    x_values[i] = x
    y_values[i] = y

#Printing results
print('The longest flight is achieved for v0=', u_set[x_max_index[0]][0,0], 'and theta0=', round(u_set[x_max_index[0]][0,1], 3), 'at', round(x_max, 3), 'metres.')

#plotting the trajetory with the largest x value
lab = 'v0='+ str(u_set[x_max_index[0]][0,0]) + ', theta0=' + str(round(u_set[x_max_index[0]][0,1], 3))
fig, ax=plt.subplots()
plt.plot(np.trim_zeros(x_values[x_max_index[0]], trim='b'), np.trim_zeros(y_values[x_max_index[0]], trim='b'), label = lab)
plt.legend()
#plt.grid(True, linetype = '--', linewidth=0.2)

print('Time taken to run this model =', time.clock() - start_time)
