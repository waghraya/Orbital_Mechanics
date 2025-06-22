import numpy as np

"""
    Inputs:
        -State OE semi_major_axis,eccentricity,true_anomaly,ascending_nodes,arg_periapsis,inclination
        -grav_param
    Outputs:
        -pos_pqw
        -vel_pqw
"""

def OE2RV(orbital_elements, grav_param):
    semi_major_axis         = orbital_elements['semi_major_axis']
    eccentricity            = orbital_elements['eccentricity']
    true_anomaly            = orbital_elements['true_anomaly']
    ascending_nodes         = orbital_elements['ascending_nodes']
    arg_periapsis           = orbital_elements['arg_periapsis']
    inclination             = orbital_elements['inclination']
    
    semi_latus_rectum       = semi_major_axis*(1-np.POW(eccentricity,2))
    # Defining position
    pos_pqw                 = np.array([
                        [semi_latus_rectum*np.cos(true_anomaly)/(1+eccentricity*np.cos(true_anomaly))],
                        [semi_latus_rectum*np.sin(true_anomaly)/(1+eccentricity*np.cos(true_anomaly))],
                        [0]
        ])
    # pos edge cases 
        # Circular equatorial (RAAN = 0, arg_peri = 0, f = lambda true?)
        # Circular inclined (arg_peri = 0 and f = u?)
    # Defining velocity
    vel_pqw                 = np.array([
                        [-np.sqrt(grav_param/semi_latus_rectum)*np.sin(true_anomaly)],
                        [np.sqrt(grav_param/semi_latus_rectum)*(eccentricity+np.cos(true_anomaly))],
                        [0]
        ])
    # vel edge cases
        # elliptical equatorial (arg_peri = arg_peri_true)