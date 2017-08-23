# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 07:59:10 2017

@author: villegas
"""

import numpy as np
import matplotlib.pyplot as plt

T = 100.
dt = 0.1
N = int(T/dt) + 1
t_space = np.linspace(0., T, N)

ms = 50.
g = 9.81
rho = 1.091
A = np.pi * 0.5**2
Cd = 0.15
mp0 = 100.
ve = 325.

#Definition of the burning rate. Truncated function
mdot = np.zeros(N)
mdot[:int(5./dt)] = 20.

#Initial conditions
h0 = 0.
v0 = 0.
u  = np.array([h0, v0])
mp = mp0
mc = 0.

def f(u):
    return np.array([u[1], -g + (mdot[n]*ve - 0.5*rho*u[1]*abs(u[1])*A*Cd)/(ms + mp)])

def euler_step(u, f, dt):
    return u + dt * f(u)

sol_set = sol_set = np.zeros((N,3))
sol_set[0][0] = t_space[0]
sol_set[0][1] = h0
sol_set[0][2] = v0


for n, t in enumerate(np.delete(t_space, -1)):
    mp = mp0  - np.sum(mdot[:n]*dt)
    
    u = euler_step(u, f, dt)
    
    sol_set[n+1][0] = t+dt
    sol_set[n+1][1] = u[0]
    sol_set[n+1][2] = u[1]
    
    if u[0] < 0:
        break
    

#for i in range(N-1):
    #mp = mp0 - np.sum(mdot[:i]*dt)
    #u = u + dt * np.array([u[1], -g + (mdot[i]*ve - 0.5*rho*u[1]*abs(u[1])*A*Cd)/(ms + mp)])
    #h[i+1] = u[0]
    #v[i+1] = u[1]