import numpy as np
import configparser
from utils.constants import omega_Earth, grav_param_Earth, radius_Earth
from utils.Rotations.SEZ2ECI import SEZ2ECI as S2E
from .base import PreliminaryOD

"""
Inputs:
    phi = 
    lam = 
    rho = scalar distance from ground station to satellite
    z = altitude of ground station relative to sea level (km)
    sig = 
    beta = 
    rho_dot = rate of change of rho with respect to time (km/s)
    sig_dot = rate change of sig with respect to time (deg/s)
    beta_dot = rate change of beta with respect to time (deg/s)
"""

class RangeAngleOD(PreliminaryOD):
    def __init__(self, sat_name, sat_loc, ground_station_loc, sim_type):
        self.sim_type           = sim_type
        self.sat_name           = sat_name
        super().__init__(sat_name, sim_type)
        self.phi                = ground_station_loc['phi']
        self.lam                = ground_station_loc['lam']
        self.altitude           = ground_station_loc['altitude']
        self.sat_range          = sat_loc['sat_range']
        self.sig                = sat_loc['sig']
        self.beta               = sat_loc['beta']
        self.sat_range_rate     = sat_loc['sat_range_rate']
        self.sig_rate           = sat_loc['sig_rate']
        self.beta_rate          = sat_loc['beta_rate']
    
    @classmethod
    def import_config(cls, config_path='config/RangeAngleOD_config.ini'):
        config = configparser.ConfigParser()
        config.read(config_path)
        
        sat_loc                                 = {}
        ground_station_loc                      = {}
        
            # ---------take inputs-----------
        sat_name                                = config['Satellite']['sat_name']
        sim_type                                = config['Scenario']['sim_type']
        ground_station_loc['phi']               = config["Ground Station"]["phi"]
        ground_station_loc['lam']               = config["Ground Station"]["lam"]
        ground_station_loc['altitude']          = config["Ground Station"]["altitude"]
        
        sat_loc['sat_range']                    = config["Satellite"]["sat_range"]
        sat_loc['sig']                          = config["Satellite"]["sig"]
        sat_loc['beta']                         = config["Satellite"]["beta"]
        sat_loc['sat_range_rate']               = config["Satellite"]["sat_range_rate"]
        sat_loc['sig_rate']                     = config["Satellite"]["sig_rate"]
        sat_loc['beta_rate']                    = config["Satellite"]["beta_rate"]

        return cls(sat_name,sat_loc,ground_station_loc,sim_type)
    
    
    def run(self):
        print("Running RangeAngleOD method")
        # Implementation of RangeAngleOD
        
    #---------------------[finding r1,v1]---------------------------
        sat_range_SEZ           = self.sat_range*np.array([
                                    [-np.cos(self.sig)*np.cos(self.beta)],
                                    [np.cos(self.sig)*np.sin(self.beta)],
                                    [np.sin(self.sig)]
                                    ])
        sat_range_ECI           = S2E(sat_range_SEZ,self.lam,self.phi)
        sat_range_rate_SEZ_SEZ  = self.sat_range_rate*np.array([
                                    [-np.cos(self.sig)*np.cos(self.beta)],
                                    [np.cos(self.sig)*np.sin(self.beta)],
                                    [np.sin(self.sig)]
                                    ]) + self.sat_range*np.array([
                                        [(self.sig_rate*(np.sin(self.sig)*np.cos(self.beta)))+(self.beta_rate*np.cos(self.sig)*np.sin(self.beta))],
                                        [(-self.sig_rate*np.sin(self.sig)*(np.sin(self.beta)))+(self.beta_rate*np.cos(self.sig)*np.cos(self.beta))],
                                        [(self.sig_rate*np.cos(self.sig))]
                                    ])
        sat_range_rate_SEZ_ECI  = S2E(sat_range_rate_SEZ_SEZ,self.lam,self.phi)
        site_pos_SEZ            = np.array([[0],[0],[radius_Earth]])
        site_pos_ECI            = S2E(site_pos_SEZ,self.lam,self.phi)
        pos1_ECI                = site_pos_ECI + sat_range_ECI
        vel1_ECI                = sat_range_rate_SEZ_ECI + np.cross(omega_Earth.T,pos1_ECI.T).T

    #---------------------[r1,v1 -> OE1]----------------------------
        orbitalElements = self.RV2OE(pos1_ECI, vel1_ECI, grav_param_Earth)
        eccentricity = orbitalElements['Eccentricity']
        semi_major_axis = orbitalElements['Semi Major Axis']
        true_anomaly = orbitalElements['True Anomaly']
        RAAN = orbitalElements['RAAN']
        arg_periapsis = orbitalElements['Argument of Periapsis']
        inclination = orbitalElements['Inclination']

        E_anomaly = 2*np.atan(np.tan(true_anomaly/2)*np.sqrt((1-np.linalg.norm(eccentricity))/(1+np.linalg.norm(eccentricity))))
        # Kepler's Equation
        M_anomaly = E_anomaly - np.linalg.norm(eccentricity)*np.sin(E_anomaly)
        mean_motion = np.sqrt(grav_param_Earth/semi_major_axis**3)
        orbitalPeriod = 2*np.pi/mean_motion
        return orbitalElements, orbitalPeriod