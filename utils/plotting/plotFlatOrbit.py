# -*- coding: utf-8 -*-
"""
Created on Sun Jul 13 12:16:25 2025

@author: waghr


Dev notes:
    Plotting orbit based off given OrbitalElements in the PQW frame
"""
import matplotlib.pyplot as mplot
import numpy as np

def plotFlatOrbit(orbitalElements):
    theta = np.linspace(0,2*np.pi,500)
    semi_major_axis = orbitalElements['Semi Major Axis']
    eccentricity = orbitalElements['Eccentricity']
    true_anomaly = orbitalElements(['True Anomaly'])
    current_radial_pos = semi_major_axis * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos(true_anomaly))
    current_x_pos = current_radial_pos*np.cos(true_anomaly)
    current_y_pos = current_radial_pos*np.sin(true_anomaly)
    radial_pos = (semi_major_axis * (1 - eccentricity ** 2)) / (1 + eccentricity * np.cos(theta))

    # convert from radial to cartesian
    x_pos = radial_pos * np.cos(theta)
    y_pos = radial_pos * np.sin(theta)

    # Plotting
    mplot.figure(figsize=(6, 6))
    mplot.plot(x_pos, y_pos, label='Orbit Path')
    mplot.plot([0], [0], 'bo', label='Earth')
    mplot.plot([semi_major_axis*(1-eccentricity)], [0], 'ro', label='Periapsis')
    mplot.plot([-semi_major_axis*(1+eccentricity)], [0], 'ro', label='Apoapsis')
    mplot.plot([current_x_pos], [current_y_pos], 'ko', label='Current Position')
    
    # Plot settings
    mplot.axis('equal')
    mplot.title('2D Orbital Plot (Orbital Plane)')
    mplot.xlabel('x (km)')
    mplot.ylabel('y (km)')
    mplot.grid(True)
    mplot.legend()
    mplot.show()
OE = {'Semi Major Axis': 7000, 'Eccentricity': 0.1}
plotFlatOrbit(OE)