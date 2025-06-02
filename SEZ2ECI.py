import numpy as np

"""
Converts vectors from SEZ (South, East, Zenith) coordinate frame to the ECI (Earth-Centered Inertial) coordinate frame
Inputs: Initial SEZ vector, lam (), phi () 
Outputs: Converted ECI vector

"""

def SEZ2ECI(vec_SEZ, lam, phi):
    rot1 = np.array([
                    [np.cos(lam), -np.sin(lam), 0], 
                    [np.sin(lam), np.cos(lam), 0], 
                    [0,0,1]
                    ])
    rot2 = np.array([
                    [np.sin(phi), 0, np.cos(phi)],
                    [0,1,0],
                    [-np.cos(phi),0,np.sin(phi)]
                    ])
    vec_ECI = rot1@rot2@vec_SEZ
    return(vec_ECI)