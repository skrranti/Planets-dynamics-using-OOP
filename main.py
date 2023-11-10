# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:11:14 2023

@author: santi
"""
import numpy as np
from SSmodule import astro
import SSmodule as ss

#%% Define the time array and the planets objects from the class astro

t_ini = 0 #ms
t_fin = 1000#ms
dt = 1e-2 #ms
t = np.arange(t_ini,t_fin+dt,dt)

mercury = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
venus = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
earth = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
mars = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
jupyter = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
saturn = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
uranus = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
neptune = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
pluto = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)
sun = astro(np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),np.zeros(len(t)),0)

solarsystem = [mercury,venus,earth,mars,jupyter,saturn,uranus,neptune,pluto,sun] #creat the list of planets


#give initial conditions
i = 0
for astros in solarsystem: #initial conditions
    astros.x[0] = ss.d_[i]
    astros.y[0] = 0
    astros.vx[0] = 0
    astros.vy[0] = ss.ov_[i]
    astros.m = ss.m_[i]
    i = i+1
    
#calculate initial accelerations given initial positions of the planets
for astros in solarsystem: #calculate initial accelerations thanks to initial positions
    astros.ginteraction(solarsystem, 0)
   

for i in range(len(t)-1): #apply verlet algorithm to evolve 
    print(i)
    for astros in solarsystem:
        astros.taylorposition(dt, i)
    
    for astros in solarsystem:
        astros.ginteraction(solarsystem, i+1)
        astros.taylorvelocity(dt, i)
    
    

#%% PLOT animation

import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

fig, ax =plt.subplots()
ax.axis([-5,5,-5,5]) ## change the axis limits to a closer view of the nearest planets orbits
ax.set_aspect("equal")


mercury_ax, = ax.plot(0,1,color='gray', marker='o', label='mercury')
venus_ax, = ax.plot(0,1,color='black', marker='o', label=' venus')
earth_ax, = ax.plot(0,1,color='blue', marker='o', label='earth')
mars_ax, = ax.plot(0,1,color='red', marker='o', label='mars')
jupyter_ax, = ax.plot(0,1,color='brown', marker='o', label='jupyter')
saturn_ax, = ax.plot(0,1,color='purple', marker='o', label='saturn')
uranus_ax, = ax.plot(0,1,color='green', marker='o', label='uranus')
neptune_ax, = ax.plot(0,1,color='crimson', marker='o', label='neptune')
pluto_ax, = ax.plot(0,1,color='springgreen', marker='o', label='pluto')
sun_ax, = ax.plot(0,1,color='orange', marker='o', label='sun')

axiss = [mercury_ax, venus_ax, earth_ax, mars_ax, jupyter_ax, saturn_ax, uranus_ax, neptune_ax, pluto_ax, sun_ax]

ax.legend()
def animate(i):
    for astros, axis in zip(solarsystem,axiss):
        axis.set_data(astros.x[i], astros.y[i])
    
   
    

ani = FuncAnimation(fig, animate, frames=len(t), interval=20, repeat=False)
plt.show()







        
        

