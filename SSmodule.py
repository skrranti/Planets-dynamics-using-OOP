# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:51:36 2023

@author: santi
"""

import numpy as np 

def distancia(xi,yi,xj,yj): #create a function to calculate the euclidean distance in a plane given coordinates x and  y of two points
    f = np.sqrt( (xi-xj)**2 + (yi-yj)**2 )
    return f*f*f

class astro: # we create the class astro 
    def __init__(self, x_array, y_array, vx_array, vy_array, ax_array, ay_array, m): # the atributes are the position , velocity and acceleration coordinates of each planet
        self.x = x_array                                                             
        self.y = y_array
        self.vx = vx_array
        self.vy = vy_array
        self.ax = ax_array
        self.ay = ay_array
        self.m = m
        
    def taylorposition(self,h,i): #this method calculates the position of a planet at a time instant given the position , velocity and acceleration at the previous time instant
        self.x[i+1] = self.x[i] + h*self.vx[i] + h*h*0.5*self.ax[i]
        self.y[i+1] = self.y[i] + h*self.vy[i] + h*h*0.5*self.ay[i]
        
    def ginteraction(self,system,i): ## calculates the acceleration at a time instant of a planet interacting with a group of planets from the class astro whose positions are known
        fx = 0
        fy = 0
        for astros in system:
            if (astros!=self):
                addx = -astros.m*(self.x[i]-astros.x[i])/(distancia(self.x[i],self.y[i],astros.x[i],astros.y[i]))
                addy = -astros.m*(self.y[i]-astros.y[i])/(distancia(self.x[i],self.y[i],astros.x[i],astros.y[i]))
                fx = fx + addx
                fy = fy + addy
        self.ax[i] = fx
        self.ay[i] = fy
        
    def taylorvelocity(self,h,i): # this method calculates the velocity of a planet given the position and acceleration of the planet in the previous time instant but also the acceleration in the same time instant
        self.vx[i+1] = self.vx[i] + 0.5*h*(self.ax[i]+self.ax[i+1])
        self.vy[i+1] = self.vy[i] + 0.5*h*(self.ay[i]+self.ay[i+1])


#%% PARAMETERS


Ms = 1.99e30 ## Kg
G = 6.67e-11 ## Nm^2/Kg^2

c = 1.496e11 ## m

t_rfactor = np.sqrt(G*Ms/c**3) 

## mass of each planet 

m = np.zeros(10)

m[0] = 0.330   #mercury
m[1] = 4.87    #venus
m[2] = 5.97    #earth
m[3] = 0.642   #mars
m[4] = 1899    #jupyter
m[5] = 568     #saturn
m[6] = 86.8    #uranus
m[7] = 102     #neptune 
m[8] = 0.0125  #pluto
m[9] = Ms*1e-24 #sun

m = m*1e24

m_ = m/Ms ## reescalated mass

## distance to the sun, we will take them as initial conditions por the coordinate x

d = np.zeros(10)

d[0] = 57.9
d[1] = 108.2
d[2] = 149.6
d[3] = 227.9
d[4] = 778.6
d[5] = 1433.5
d[6] = 2872.5
d[7] = 4495.1
d[8] = 5870.0
d[9] = 0

d = d*1e9

d_ = d/c ## distance reescalated

## orbital velocities that will be taken as initial conditions for coordinate y of velocity

ov = np.zeros(10)

ov[0] = 47.9
ov[1] = 35.0
ov[2] = 29.8
ov[3] = 24.1
ov[4] = 13.1
ov[5] = 9.7
ov[6] = 6.8
ov[7] = 5.4
ov[8] = 4.7
ov[9] = 0

ov = ov*1000

ov_ = ov/(c*t_rfactor) ## reescalated velocity

