import numpy as np
import pandas as pd
from utils.constants import grav_param_Earth


class Propagator:
    def __init__(self, orbitalPeriod, propagator_type, initial_pos_ECI, initial_vel_ECI):
        self.propagator_type = 'N/A'
        self.orbitalPeriod = orbitalPeriod
        self.state_vector_ECI = np.hstack(initial_pos_ECI, initial_vel_ECI)
        self.propagatedStates = pd.DataFrame(columns=[
            'time',
            'position_I', 'position_J', 'position_K',
            'velocity_I', 'velocity_J', 'velocity_K',
            'acceleration_I', 'acceleration_J', 'acceleration_K',
            'semi_major_axis', 'eccentricity', 'true_anomaly', 'inclination', 'raan', 'arg_periapsis'])
        self.state = np.array([[initial_pos_ECI],
                               [initial_vel_ECI],
                               [-grav_param_Earth/(np.linalg.norm(initial_pos_ECI))*initial_pos_ECI]])
    def Propagate(self):
        print('No propagator selected')

    @staticmethod
    def EulerSolver(initial_pos_km, initial_vel_kms, initial_acc_kms2, timeStep, grav_param):
        r_norm = np.linalg.norm(initial_pos_km)
        acc_kms2 = -grav_param * initial_pos_km/ (r_norm ** 3)
        vel_kms_ECI = initial_vel_kms + acc_kms2 * timeStep
        pos_km_ECI = initial_pos_km + initial_vel_kms * timeStep
        return [pos_km_ECI, vel_kms_ECI, acc_kms2]

    @staticmethod
    def two_body_ODE(self):
        acc_kms2 = -grav_param_Earth/(np.linalg.norm(self.stat()))*pos_km
        return [vel_kms, acc_kms2]
