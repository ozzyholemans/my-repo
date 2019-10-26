#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:11:10 2019

@author: oscarholemans
"""

import numpy as np
import matplotlib.pyplot as plt

#define same constants as before 
G=6.67e-11
earth_mass=5.972e24
moon_mass=7.34767309e22
earth_radius=6.371e6
moon_radius=1.737e6

#new acceleration due to two moving bodies
def accel_1(r1,r2):
    A1 = -(earth_mass*moon_mass*G*(r1-r2))/np.linalg.norm(r1-r2)**3
    return A1

def accel_2(r1,r2):
    A2 = -(earth_mass*moon_mass*G*(r2-r1))/np.linalg.norm(r2-r1)**3
    return A2
    

#new function includes both positional arguments 
def runge_kutta(function,position1,position2,velocity,dt):
    
    K1x=velocity
    K1v=function(position1,position2)
    print(K1x)
    
    K2x=velocity+(dt*K1v/2)
    K2v=function(position1+dt*K1x/2,position2+dt*K1x/2)
    print(K2x)
    
    K3x=velocity+(dt*K2v/2)
    K3v=function(position1+dt*K2x/2,position2+dt*K2x/2)
    print(K3x)
    
    K4x=velocity+(dt*K3v)
    K4v=function(position1+dt*K3x,position2+dt*K3x/2)
    print(K4x)
    
    Pos_step=dt*(K1x+2*K2x+2*K3x+K4x)/6
    Vel_step=dt*(K1v+2*K2v+2*K3v+K4v)/6
    
    return [Pos_step,Vel_step]

MyInput = '0'
while MyInput != 'q':
    MyInput = input('Enter a choice, "2 body system (Earth and moon) (a)","3 body system (Earth and moon with the Sun) (b)" or "(q)" to quit: ')
    
    #set empty arrays for values to be assigned to
    time_vals = []
    x1_vals = []
    y1_vals = []
    x2_vals = []
    y2_vals = []
    
    if MyInput == 'a':
        print('You have chosen part (a)')
        runtime=float(input('Please now choose the runtime of the simulation in seconds (float):'))
        
        time = 0
        dt = 1
        #initialise problem with initial location of earth and moon, and velocity
        earth_position = np.array([0.0,0.1])
        earth_velocity = np.array([500.0,100.0])
        moon_position = np.array([0.0,384400e3])
        moon_velocity = np.array([-500.0,-100.0])
      
        while time<runtime:
            x1_vals += [earth_position.item(0)]
            y1_vals += [earth_position.item(1)]
            x2_vals += [moon_position.item(0)]
            y2_vals += [moon_position.item(1)]
    
            time_vals += [time]
            time += dt
            #perform the runge-kutta method for the earth 
            position_velocity1 = runge_kutta(accel_1,earth_position,moon_position,earth_velocity,dt)
            earth_position += position_velocity1[0]
            earth_velocity += position_velocity1[1]
            #now perform runge-kutta method for the moon
            position_velocity2 = runge_kutta(accel_2,moon_position,earth_position,moon_velocity,dt)
            moon_position += position_velocity2[0]
            moon_velocity += position_velocity2[1]
            #print values so can track where the function is blowing up
            print(earth_position)
            print(moon_position)
        
        #plot path of earth and moon 
        plt.title('Plot of moon and the Earth')
        plt.xlabel('x(m)')
        plt.ylabel('y(m)')
        plt.plot(x1_vals,y1_vals,'r',label='Path of earth')
        plt.plot(x2_vals,y2_vals,'c',label='Path of moon')
        plt.legend()
        plt.show()
        plt.clf()
        
    #this option was to troubleshoot my code and see what was going wrong
    elif MyInput == 'b':
        time = 0
        dt = 10
        earth_position = np.array([0.0,0.1])
        earth_velocity = np.array([0.0,1.0])
        moon_position = np.array([0.0,384400e3])
        moon_velocity = np.array([1.0,0.0])    
        runge_kutta(accel_1,earth_position,moon_position,earth_velocity,dt)
        runge_kutta(accel_2,earth_position,moon_position,earth_velocity,dt)
        print(accel_1(earth_position,moon_position))
            