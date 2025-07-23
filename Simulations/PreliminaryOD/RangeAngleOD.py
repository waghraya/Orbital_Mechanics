import numpy as np
import math
import configparser
from utils.plotting.plotFlatOrbit import plotFlatOrbit
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
    def __init__(self, sat_name, sat_loc, ground_station_loc, sim_type, postProcessingDict):
        self.sim_type           = sim_type
        self.sat_name           = sat_name
        super().__init__(sat_name, sim_type)
        self.phi                = float(np.pi*float(ground_station_loc['phi'])/180)
        self.lam                = float(np.pi*float(ground_station_loc['lam'])/180)
        self.altitude           = float(ground_station_loc['altitude'])
        self.sat_range          = float(sat_loc['sat_range'])
        self.sigma              = float(np.pi*float(sat_loc['sigma'])/180)
        self.beta               = float(np.pi*float(sat_loc['beta'])/180)
        self.sat_range_rate     = float(sat_loc['sat_range_rate'])
        self.sigma_rate         = float(np.pi*float(sat_loc['sigma_rate'])/180)
        self.beta_rate          = float(np.pi*float(sat_loc['beta_rate'])/180)
        self.bool_plotFlatOrbit = bool(postProcessingDict['bool_plotFlatOrbit'])
    
    @classmethod
    def import_config(cls, config_path='config/RangeAngleOD_config.ini'):
        config = configparser.ConfigParser()
        config.read(config_path)
        
        sat_loc                                     = {}
        ground_station_loc                          = {}
        postProcessingDict                          = {}
            # ---------take inputs-----------
        sat_name                                    = config['Satellite']['sat_name']
        sim_type                                    = config['Scenario']['sim_type']
        postProcessingDict['bool_plotFlatOrbit']    = config['Scenario']['bool_plotFlatOrbit']
        ground_station_loc['phi']                   = config["Ground Station"]["phi"]
        ground_station_loc['lam']                   = config["Ground Station"]["lam"]
        ground_station_loc['altitude']              = config["Ground Station"]["altitude"]
        
        sat_loc['sat_range']                        = config["Satellite"]["sat_range"]
        sat_loc['sigma']                            = config["Satellite"]["sigma"]
        sat_loc['beta']                             = config["Satellite"]["beta"]
        sat_loc['sat_range_rate']                   = config["Satellite"]["sat_range_rate"]
        sat_loc['sigma_rate']                       = config["Satellite"]["sigma_rate"]
        sat_loc['beta_rate']                        = config["Satellite"]["beta_rate"]

        return cls(sat_name,sat_loc,ground_station_loc,sim_type,postProcessingDict)
    
    
    def run(self):
        print("Running RangeAngleOD method")
        # Implementation of RangeAngleOD
    #---------------------[finding r1,v1]---------------------------
        sat_range_SEZ           = self.sat_range*np.array([
                                    [-np.cos(self.sigma)*np.cos(self.beta)],
                                    [np.cos(self.sigma)*np.sin(self.beta)],
                                    [np.sin(self.sigma)]
                                    ])
        sat_range_ECI           = S2E(sat_range_SEZ,self.lam,self.phi)
        sat_range_rate_SEZ_SEZ  = self.sat_range_rate*np.array([
                                    [-np.cos(self.sigma)*np.cos(self.beta)],
                                    [np.cos(self.sigma)*np.sin(self.beta)],
                                    [np.sin(self.sigma)]
                                    ]) + self.sat_range*np.array([
                                        [(self.sigma_rate*(np.sin(self.sigma)*np.cos(self.beta)))+(self.beta_rate*np.cos(self.sigma)*np.sin(self.beta))],
                                        [(-self.sigma_rate*np.sin(self.sigma)*(np.sin(self.beta)))+(self.beta_rate*np.cos(self.sigma)*np.cos(self.beta))],
                                        [(self.sigma_rate*np.cos(self.sigma))]
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
        e = np.linalg.norm(eccentricity)
        E_anomaly = 2 * math.atan(math.tan(true_anomaly / 2) * math.sqrt((1 - e) / (1 + e)))
        # Kepler's Equation
        M_anomaly = E_anomaly - np.linalg.norm(eccentricity)*np.sin(E_anomaly)
        mean_motion = np.sqrt(grav_param_Earth/semi_major_axis**3)
        orbitalPeriod = 2*np.pi/mean_motion

        # Post-processing stuff
        if self.bool_plotFlatOrbit:
            plotFlatOrbit(orbitalElements)

        return orbitalElements, orbitalPeriod