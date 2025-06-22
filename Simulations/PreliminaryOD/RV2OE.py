import numpy as np

"""
    Inputs:
        -State pos,vel,grav_param
        -In ECI frame
    Outputs:
        -State OE semi_major_axis,eccentricity,true_anomaly,ascending_nodes,arg_periapsis,inclination
"""

def RV2OE(pos,vel,grav_param):
    norm_pos                = np.linalg.norm(pos)
    norm_vel                = np.linalg.norm(vel)
    spec_ang_momentum       = np.cross(pos,vel)
    K_vec                   = np.array([
                                [0],
                                [0],
                                [1]
                                ]) # i,j,k = <0,0,1>
    line_nodes              = np.cross(K_vec,spec_ang_momentum)
    spec_energy             = 0.5*np.POW(norm_vel,2) - grav_param/norm_pos
    eccentricity            = (spec_energy*pos)-(np.dot(pos,vel)*vel)/grav_param
    
    if (np.linalg.norm(eccentricity) != 0):
        semi_major_axis     = -grav_param/spec_energy
    else:
        #Parabolic orbit
        semi_major_axis     = float('inf')
        
    inclination             = np.acos(spec_ang_momentum[2,0]/np.linalg.norm(spec_ang_momentum))
    
    ascending_nodes         = np.acos(line_nodes[0,0]/np.linalg.norm(line_nodes))
    if (line_nodes[1,0] < 0): 
        ascending_nodes     = 2*np.pi - ascending_nodes
    
    arg_periapsis           = np.acos(np.dot(line_nodes,eccentricity)/(np.linalg.norm(line_nodes)*np.linalg.norm(eccentricity)))
    if (eccentricity[2,0] < 0):
        arg_periapsis       = 2*np.pi - arg_periapsis
    
    true_anomaly            = np.acos(np.dot(eccentricity,pos)/(np.linalg.norm(eccentricity)*np.linalg.norm(pos)))
    if (np.dot(pos,vel) < 0):
        true_anomaly        = 2*np.pi - true_anomaly
    
    #---Special cases---
    #Elliptical equatorial
    #Circular inclined
    #Circular equatorial
    
    OE = [semi_major_axis,eccentricity,true_anomaly,ascending_nodes,arg_periapsis,inclination]
    return OE