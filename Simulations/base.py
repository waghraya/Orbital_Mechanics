# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 14:20:25 2025

@author: waghr
"""
import numpy as np
import math
class Simulation:
    def __init__(self, sat_name, sim_type):
        self.sim_type = sim_type
        self.sat_name = sat_name
    def run(self):
        raise NotImplementedError("Subclasses must implement this method.")

    @staticmethod
    def RV2OE(pos, vel, grav_param):
        norm_pos = np.linalg.norm(pos)
        norm_vel = np.linalg.norm(vel)
        spec_ang_momentum = np.cross(pos.T, vel.T).T
        K_vec = np.array([
            [0],
            [0],
            [1]
        ])  # i,j,k = <0,0,1>
        line_nodes = np.cross(K_vec.T, spec_ang_momentum.T).T
        spec_energy = 0.5 * norm_vel**2 - grav_param / norm_pos
        eccentricity = (spec_energy * pos) - (float(pos.T @ vel) * vel) / grav_param
        if np.linalg.norm(eccentricity) != 0:
            semi_major_axis = -grav_param / spec_energy
        else:
            # Parabolic orbit
            semi_major_axis = float('inf')
        inclination = math.acos(spec_ang_momentum[2, 0] / np.linalg.norm(spec_ang_momentum))
        ascending_nodes = math.acos(line_nodes[0, 0] / np.linalg.norm(line_nodes))
        if line_nodes[1, 0] < 0:
            ascending_nodes = 2 * np.pi - ascending_nodes

        arg_periapsis = math.acos(
            float(line_nodes.T @ eccentricity) / (np.linalg.norm(line_nodes) * np.linalg.norm(eccentricity)))
        if eccentricity[2, 0] < 0:
            arg_periapsis = 2 * np.pi - arg_periapsis

        true_anomaly = math.acos(float(eccentricity.T @ pos) / (np.linalg.norm(eccentricity) * np.linalg.norm(pos)))
        if float(pos.T @ vel) < 0:
            true_anomaly = 2 * np.pi - true_anomaly

        # ---Special cases---
        # Elliptical equatorial
        # Circular inclined
        # Circular equatorial

        orbitalElements = {'Semi Major Axis':semi_major_axis,
              'Eccentricity':eccentricity,
              'True Anomaly':true_anomaly,
              'RAAN':ascending_nodes,
              'Argument of Periapsis':arg_periapsis,
              'Inclination':inclination}
        return orbitalElements

    @staticmethod
    def OE2RV(orbitalElements, grav_param):
        semi_major_axis = orbitalElements['Semi Major Axis']
        eccentricity = orbitalElements['Eccentricity']
        true_anomaly = orbitalElements['True Anomaly']
        ascending_nodes = orbitalElements['RAAN']
        arg_periapsis = orbitalElements['Argument of Periapsis']
        inclination = orbitalElements['Inclination']

        semi_latus_rectum = semi_major_axis * (1 - eccentricity**2)
        # Defining position
        pos_pqw = np.array([
            [semi_latus_rectum * np.cos(true_anomaly) / (1 + eccentricity * np.cos(true_anomaly))],
            [semi_latus_rectum * np.sin(true_anomaly) / (1 + eccentricity * np.cos(true_anomaly))],
            [0]
        ])
        # pos edge cases
        # Circular equatorial (RAAN = 0, arg_peri = 0, f = lambda true?)
        # Circular inclined (arg_peri = 0 and f = u?)
        # Defining velocity
        vel_pqw = np.array([
            [-np.sqrt(grav_param / semi_latus_rectum) * np.sin(true_anomaly)],
            [np.sqrt(grav_param / semi_latus_rectum) * (eccentricity + np.cos(true_anomaly))],
            [0]
        ])
        # vel edge cases
        # elliptical equatorial (arg_peri = arg_peri_true)
        return [pos_pqw, vel_pqw]