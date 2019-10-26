#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 14:19:35 2019

@author: oscarholemans
"""

import numpy as np
import matplotlib.pyplot as plt

#define some useful constants for use in the equations
G=6.67e-11
earth_mass=5.972e24
moon_mass=7.34767309e22
earth_radius=6.371e6
moon_radius=1.737e6

moon_position=np.array([0,384400e3])

#calculates acceleration due to earth
def e_acceleration(r):
    A = -(earth_mass*G*r)/np.linalg.norm(r)**3
    return A

#calculates acceleration due to earth and moon
def e_m_acceleration(r):
    A_earth = -(earth_mass*G*r)/np.linalg.norm(r)**3
    A_moon = -(moon_mass*G*(r-moon_position))/(np.linalg.norm(r-moon_position)**3)
    return A_earth+A_moon

#calculates the potential energy due to either the earth or the earth and moon
def potential_energy(r, moon='False'):
    earth_potential = -(G*earth_mass*rocket_mass)/np.linalg.norm(r)
    moon_potential = -(G*moon_mass*rocket_mass)/np.linalg.norm(r-moon_position-moon_radius)
    if moon=='False':
        return earth_potential
    elif moon=='True':
        return earth_potential+moon_potential
    
#calculates the kinetic energy for a given mass and velocity
def kinetic_energy(mass,velocity):
    ke = 0.5*mass*np.linalg.norm(velocity)**2
    return ke

#calculates the escape velocity for a position from the earth
def escape_velocity(position):
    escvel = np.sqrt(2*G*earth_mass/np.linalg.norm(position))
    return escvel

#calculates the velocity for a circular orbit at a position
def circular_orbit(position):
    circvel = np.sqrt(G*earth_mass/np.linalg.norm(position)) 
    return circvel

#plots circular earth
def earth_plot():
    global earth_x
    global earth_y
    earth_x=[]
    earth_y=[]
    circle = np.linspace(0,2*np.pi,100)
    
    for i in circle:
        earth_x+=[earth_radius*np.cos(i)]
        earth_y+=[earth_radius*np.sin(i)]
        
    return

#plots circular moon
def moon_plot():
    global moon_x
    global moon_y
    moon_x=[]
    moon_y=[]
    circle = np.linspace(0,2*np.pi,100)
    
    for i in circle:
        moon_x+=[moon_radius*np.cos(i)]
        moon_y+=[384400e3+moon_radius*np.sin(i)]
        
    return

#numerically solves the function inputted for a time-step
def runge_kutta(function,position,velocity,dt):
    #first step
    K1x=velocity
    K1v=function(position)
    #second step
    K2x=velocity+(dt*K1v/2)
    K2v=function(position+dt*K1x/2)
    #third step
    K3x=velocity+(dt*K2v/2)
    K3v=function(position+dt*K2x/2)
    #fourth step
    K4x=velocity+(dt*K3v)
    K4v=function(position+dt*K3x)
    #total step
    x_step=dt*(K1x+2*K2x+2*K3x+K4x)/6
    v_step=dt*(K1v+2*K2v+2*K3v+K4v)/6
    #returns an array with the values in
    return [x_step,v_step]


#menu to allow user to choose which part of the program they want to run
MyInput = '0'
while MyInput != 'q':
    MyInput = input('Enter a choice, "Rocket trajectory around the Earth (a)","Rocket mission to take photos of the moon (b)" or "(q)" to quit: ')
    #define empty arrays to assign answers to
    time_vals = []
    x_vals = []
    y_vals = []
    kinetic_vals = []
    potential_vals = []
    total_energy_vals = []
    moon_distance = []
    if MyInput == 'a':
        print('You have chosen part (a)')
        #allow user to choose mass of rocket, initial height, initial velocity,runtime etc
        rocket_mass=float(input('Please choose the mass of the rocket in kg (float):'))
        height=float(input("Please choose the starting height of the rocket above the Earth's surface in metres (float):"))
        initial_position = height + earth_radius
        rocket_position = np.array([0,initial_position])
        
        print('From this position, the escape velocity is: ' + str(escape_velocity(initial_position)))
        print('The velocity required to achieve a circular orbit is: ' + str(circular_orbit(initial_position)))
        initial_velocity=float(input('Please now choose what the initial velocity of the rocket will be in m/s (float):'))
        rocket_velocity = np.array([initial_velocity,0])
        
        dt = 1
        time = 0
        runtime=float(input('Please now choose the runtime of the simulation in seconds (float):'))
        
        #loop stops if time exceeds runtime or the rocket 'crashes' into earth or moon
        while time<runtime and np.linalg.norm(rocket_position)>=earth_radius:
            time_vals += [time]
            #increase time step by 1
            time += dt
            #calcluate the velocity and position for the time step with the runge-kutta method 
            position_velocity = runge_kutta(e_acceleration,rocket_position,rocket_velocity,dt)
            #update rocet position and velocity arrays with the new position as calculated
            rocket_position += position_velocity[0]
            rocket_velocity += position_velocity[1]
            #update x and y values with the values from the rocket position array
            x_vals += [rocket_position.item(0)]
            y_vals += [rocket_position.item(1)]
            #calculate the normal distance so can check whether the rocket has crashed
            norm_distance = np.linalg.norm(rocket_position)
            #now calculate the kinetic and potential energy at current position and velocity
            kinetic_vals += [kinetic_energy(rocket_mass,rocket_velocity)]
            potential_vals += [potential_energy(norm_distance)]
            #therefore calculate the total energy 
            total_energy_vals += [kinetic_energy(rocket_mass,rocket_velocity) + potential_energy(norm_distance)]
        
        #plot the earth
        earth_plot()
        
        #plot the path of the rocket
        plt.title('Path of rocket around the Earth')
        plt.xlabel('x(m)')
        plt.ylabel('y(m)')
        plt.plot(x_vals,y_vals,'r',label='Path of rocket')
        plt.plot(earth_x,earth_y,'g',label='Earth')
        plt.legend()
        plt.show()
        
        #plot the kinetic, potential and total energy throughout the simulation
        plt.title('Energy plots')
        plt.xlabel('time(s)')
        plt.ylabel('energy(J)')
        plt.plot(time_vals,kinetic_vals,'r',label='Kinetic energy')
        plt.plot(time_vals,potential_vals,'g',label='Potential energy')
        plt.plot(time_vals,total_energy_vals,'b',label='Total energy')
        plt.legend()
        plt.show()
    
    elif MyInput == 'b':
        print('You have chosen part (b)')
        rocket_mass=float(input('Please choose the mass of the rocket in kg (float):'))
    
        rocket_position = np.array([0,earth_radius+6.9e6])
        rocket_velocity = np.array([115,escape_velocity(rocket_position)-135])
         
        dt = 10
        time = 0
        
        #loop will stop if the rocket crahses into either the earth or the moon, or the runtime is exceeded
        while time<700000 and np.linalg.norm(rocket_position)>=earth_radius and np.linalg.norm(rocket_position-moon_position)>=moon_radius:
            time_vals += [time]
            time += dt
            #now use the runge-kutta method but with the acceleration due to the earth and moon this time
            position_velocity = runge_kutta(e_m_acceleration,rocket_position,rocket_velocity,dt)
            #repeat same process as before to assign the values to arrays to plot
            rocket_position += position_velocity[0]
            rocket_velocity += position_velocity[1]
            
            x_vals += [rocket_position.item(0)]
            y_vals += [rocket_position.item(1)]
            
            norm_distance = np.linalg.norm(rocket_position)
            #calculates the distance to the moon so the user can know
            moon_distance += [np.linalg.norm(rocket_position-moon_position)-moon_radius]
            
            kinetic_vals += [kinetic_energy(rocket_mass,rocket_velocity)]
            potential_vals += [potential_energy(norm_distance)]
            total_energy_vals += [kinetic_energy(rocket_mass,rocket_velocity) + potential_energy(norm_distance)]
        
        #plot the earth and moon on the graph
        earth_plot()
        moon_plot()
        
        #plot the path of the rocket
        plt.title('Plot of rocket trip around the moon')
        plt.xlabel('x(m)')
        plt.ylabel('y(m)')
        plt.plot(x_vals,y_vals,'r',label='Path of rocket')
        plt.plot(earth_x,earth_y,'g',label='Earth')
        plt.plot(moon_x,moon_y,'b',label='Moon')
        plt.legend()
        plt.show()
        
        #plot the kinetic, potential and total energy throughout the simulation
        plt.title('Energy plots')
        plt.xlabel('time(s)')
        plt.ylabel('energy(J)')
        plt.plot(time_vals,kinetic_vals,'r',label='Kinetic energy')
        plt.plot(time_vals,potential_vals,'g',label='Potential energy')
        plt.plot(time_vals,total_energy_vals,'b',label='Total energy')
        plt.legend()
        plt.show()
        #tell user how close the rocket was to the moon's surface
        print('Closest approach to the moon is: ' + str(min(moon_distance)) + ' metres.')
    
    elif MyInput != 'q':
        print('This is not a valid choice')
    
print('You have chosen to finish - goodbye.')
