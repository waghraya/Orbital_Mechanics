import numpy as np
import ast
import configparser
from utils.constants import grav_param_Earth
from .base import PreliminaryOD
class GibbsOD(PreliminaryOD):
    def __init__(self,sat_name,sim_type,pos1_ECI,pos2_ECI,pos3_ECI):
        self.sat_name       = sat_name
        self.sim_type       = sim_type
        super().__init__(sat_name,sim_type)
        self.pos1_ECI       = pos1_ECI
        self.pos2_ECI       = pos2_ECI
        self.pos3_ECI       = pos3_ECI

    @classmethod
    def import_config(cls, config_path='config/GibbsOD_config.ini'):
        config = configparser.ConfigParser()
        config.read(config_path)

        sim_type         = config['Scenario']['sim_type']
        sat_name         = config['Satellite']['sat_name']
        pos1_ECI         = np.array(ast.literal_eval(config['Satellite']['pos1_ECI'])).T
        pos2_ECI         = np.array(ast.literal_eval(config['Satellite']['pos2_ECI'])).T
        pos3_ECI         = np.array(ast.literal_eval(config['Satellite']['pos3_ECI'])).T

        return cls(sat_name,sim_type,pos1_ECI,pos2_ECI,pos3_ECI)

    def run(self):
        print("Running GibbsOD method")

        norm_pos1       = np.linalg.norm(self.pos1_ECI)
        norm_pos2       = np.linalg.norm(self.pos2_ECI)
        norm_pos3       = np.linalg.norm(self.pos3_ECI)

        Z_12            = np.cross(self.pos1_ECI,self.pos2_ECI)
        Z_23            = np.cross(self.pos2_ECI,self.pos3_ECI)
        Z_31            = np.cross(self.pos3_ECI,self.pos1_ECI)

        N               = norm_pos1*Z_23+norm_pos2*Z_31+norm_pos3*Z_12
        D               = Z_12+Z_23+Z_31
        S               = (norm_pos2-norm_pos3)*self.pos1_ECI+(norm_pos3-norm_pos1)*self.pos2_ECI+(norm_pos1-norm_pos2)*self.pos3_ECI
        B               = np.cross(D,self.pos2_ECI)
        Lg              = np.sqrt(grav_param_Earth/(np.linalg.norm(N)*np.linalg.norm(D)))

        vel1_ECI        = (Lg/norm_pos1)*B + Lg*S
        vel2_ECI        = (Lg/norm_pos2)*B + Lg*S
        vel3_ECI        = (Lg/norm_pos3)*B + Lg*S

        orbitalElements = self.RV2OE(self.pos1_ECI, vel1_ECI,grav_param_Earth)
        orbitalPeriod = 2*np.pi*np.sqrt(orbitalElements['Semi Major Axis']**3/grav_param_Earth)

        return orbitalElements, orbitalPeriod