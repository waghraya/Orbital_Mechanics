import numpy as np
import configparser
from utils.constants import omega_Earth, grav_param_Earth, radius_Earth
from utils.Newton_Raphson import Newton_Raphson as NR
from utils.Rotations.SEZ2ECI import SEZ2ECI as S2E
from utils.Rotations.PQW2ECI import PQW2ECI as P2E
from utils.Rotations.PQW2ECI import ECI2PQW as E2P
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
        orbitalElements_1 = self.RV2OE(pos1_ECI, vel1_ECI, grav_param_Earth)
        eccentricity = orbitalElements_1['Eccentricity']
        semi_major_axis = orbitalElements_1['Semi Major Axis']
        true_anomaly_1 = orbitalElements_1['True Anomaly']
        RAAN_1 = orbitalElements_1['RAAN']
        arg_periapsis_1 = orbitalElements_1['Argument of Periapsis']
        inclination = orbitalElements_1['Inclination']

        E_anomaly_1 = 2*np.atan(np.tan(true_anomaly_1/2)*np.sqrt((1-np.linalg.norm(eccentricity))/(1+np.linalg.norm(eccentricity))))
        # Kepler's Equation
        M_anomaly_1 = E_anomaly_1 - np.linalg.norm(eccentricity)*np.sin(E_anomaly_1)
        mean_motion = np.sqrt(grav_param_Earth/semi_major_axis**3)
        eccentricity_norm = np.linalg.norm(eccentricity)
        # Convert r1,v1 to PQW
        pos1_PQW = E2P(pos1_ECI,RAAN_1,arg_periapsis_1,inclination)
        vel1_PQW = E2P(vel1_ECI,RAAN_1,arg_periapsis_1,inclination)

            #----------------------[OE1 -> OE2]-----------------------------
        # add perturbations for Omega and omega
        RAAN_2 = RAAN_1 # temporary
        arg_periapsis_2 = arg_periapsis_1 # temporary

        M_anomaly_2 = M_anomaly_1 + mean_motion*self.time_of_flight

        # Newton-Raphson method used to solve for E2
        E_anomaly_2 = NR(M_anomaly_2,eccentricity_norm)
        true_anomaly_2 = 2 * np.atan(np.tan(E_anomaly_2/2) * np.sqrt((1 + np.linalg.norm(eccentricity)) / (1 - np.linalg.norm(eccentricity))))
        orbitalElements_2 = {'Eccentricity': eccentricity,
                            'Semi Major Axis': semi_major_axis,
                            'True Anomaly': true_anomaly_2,
                            'RAAN': RAAN_2,
                            'Argument of Periapsis': arg_periapsis_2,
                            'Inclination': inclination}
            #---------------------[OE2 -> r2,v2]----------------------------
        [pos2_PQW, vel2_PQW] = self.OE2RV(orbitalElements_2,grav_param_Earth)
        # Final r and v converted from PQW to ECI frame
        pos2_ECI = P2E(pos2_PQW,RAAN_2,arg_periapsis_2,inclination)
        vel2_ECI = P2E(vel2_PQW,RAAN_2,arg_periapsis_2,inclination)
        pos_vec = np.array([pos1_PQW,pos2_PQW])
        vel_vec = np.array([vel1_PQW,vel2_PQW])

        '''
                    #--------------------[J2 Perturbation]---------------------------
                # Change of RAAN
                Ohm_dot = -( ((3*n*J2)/(2*( 1- e**2)**2)) * (rE/a)**2 * np.cos(i))
                Ohm = Ohm + Ohm_dot*tof
                # Change of argument of periapsis
                omega_dot = (3*n*J2)/(4*(1-e**2)**2)*((rE/a)**2)*(5*(np.cos(i)**2) - 1)
                omega = omega + omega_dot*tof
        '''