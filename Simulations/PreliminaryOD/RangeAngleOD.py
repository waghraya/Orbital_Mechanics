import numpy as np
import configparser
from utils.constants import omega_Earth, grav_param_Earth, radius_Earth
from utils.Newton_Raphson import Newton_Raphson as NR
from utils.SEZ2ECI import SEZ2ECI as S2E
from utils.PQW2ECI import PQW2ECI as P2E
from utils.PQW2ECI import ECI2PQW as E2P
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
    def __init__(self, sat_name, time_of_flight, sat_loc, ground_station_loc, sim_type):
        self.sim_type           = sim_type
        self.sat_name           = sat_name
        super().__init__(sat_name, sim_type)
        self.time_of_flight     = time_of_flight
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
        time_of_flight                          = config['Scenario']['time_of_flight']
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

        return cls(sat_name,time_of_flight,sat_loc,ground_station_loc,sim_type)
    
    
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
        K = np.array([[0],[0],[1]])
        spec_ang_momentum           = np.cross(pos1_ECI.T,vel1_ECI.T).T
        spec_ang_momentum_unit      = spec_ang_momentum/np.linalg.norm(spec_ang_momentum)
        pos1_ECI_unit                = pos1_ECI/np.linalg.norm(pos1_ECI)
        spec_energy                 = ((np.linalg.norm(vel1_ECI)**2)/2) - grav_param_Earth/(np.linalg.norm(pos1_ECI))
        line_nodes                  = np.cross(K.T,spec_ang_momentum_unit.T).T/(np.linalg.norm(np.cross(K.T,spec_ang_momentum_unit.T).T))
        eccentricity                = (np.cross(vel1_ECI.T,spec_ang_momentum.T).T/grav_param_Earth)-pos1_ECI_unit
        eccentricity_unit           = eccentricity/np.linalg.norm(eccentricity)
        semi_major_axis             = -grav_param_Earth/(2*spec_energy)
        inclination                 = float(np.acos(line_nodes.T@eccentricity_unit))

        # checking domain for omega and correcting
        if eccentricity_unit.T@K >= 0:
            arg_periapsis = float(np.acos(line_nodes.T@eccentricity_unit))
        else:
            arg_periapsis = float(2*np.pi - np.acos(line_nodes.T@eccentricity_unit))
        # atan2 autochecks quadrant
        RAAN: float = float(np.atan2(line_nodes[1],line_nodes[0]))
        # checking domain for f and correcting
        if pos1_ECI.T@vel1_ECI >= 0:
            true_anomaly_1 = np.acos(eccentricity_unit.T@pos1_ECI_unit)
        else:
            true_anomaly_1 = 2*np.pi - np.acos(eccentricity_unit.T@pos1_ECI_unit)
        E_anomaly_1 = 2*np.atan(np.tan(true_anomaly_1/2)*np.sqrt((1-np.linalg.norm(eccentricity))/(1+np.linalg.norm(eccentricity))))
        # Kepler's Equation
        M_anomaly_1 = E_anomaly_1 - np.linalg.norm(eccentricity)*np.sin(E_anomaly_1)
        mean_motion = np.sqrt(grav_param_Earth/semi_major_axis**3)
        eccentricity_norm = np.linalg.norm(eccentricity)
        # Convert r1,v1 to PQW
        pos1_PQW = E2P(pos1_ECI,RAAN,arg_periapsis,inclination)
        vel1_PQW = E2P(vel1_ECI,RAAN,arg_periapsis,inclination)

            #----------------------[OE1 -> OE2]-----------------------------

        M_anomaly_2 = M_anomaly_1 + mean_motion*self.time_of_flight

        # Newton-Raphson method used to solve for E2
        E_anomaly_2 = NR(M_anomaly_2,eccentricity_norm)
        true_anomaly_2 = 2 * np.atan(np.tan(E_anomaly_2/2) * np.sqrt((1 + np.linalg.norm(eccentricity)) / (1 - np.linalg.norm(eccentricity))))
            #---------------------[OE2 -> r2,v2]----------------------------
        pos2_norm = semi_major_axis*(1-eccentricity_norm**2)/(1+eccentricity_norm*np.cos(true_anomaly_2))
        # Final r and v determined in PQW frame
        pos2_PQW = pos2_norm*np.array([[float(np.cos(true_anomaly_2))],[float(np.sin(true_anomaly_2))],[0]])
        vel2_PQW = np.sqrt(grav_param_Earth/(semi_major_axis*(1-eccentricity_norm**2)))*np.array([[float(-np.sin(true_anomaly_2))],[eccentricity_norm+float(np.cos(true_anomaly_2))],[0]])
        # Final r and v converted from PQW to ECI frame
        pos2_ECI = P2E(pos2_PQW,RAAN,arg_periapsis,inclination)
        vel2_ECI = P2E(vel2_PQW,RAAN,arg_periapsis,inclination)
        pos_vec = np.array([pos1_PQW,pos2_PQW])
        vel_vec = np.array([vel1_PQW,vel2_PQW])

        '''
                    #--------------------[J2 Perurbation]---------------------------
                # Change of RAAN
                Ohm_dot = -( ((3*n*J2)/(2*( 1- e**2)**2)) * (rE/a)**2 * np.cos(i))
                Ohm = Ohm + Ohm_dot*tof
                # Change of argument of periapsis
                omega_dot = (3*n*J2)/(4*(1-e**2)**2)*((rE/a)**2)*(5*(np.cos(i)**2) - 1)
                omega = omega + omega_dot*tof
        '''