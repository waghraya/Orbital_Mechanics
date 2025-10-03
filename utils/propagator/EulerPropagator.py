import pandas as pd
import numpy as np
from base import Propagator
from utils.constants import grav_param_Earth

class EulerPropagator(Propagator):
    def Propagate(self, orbitalPeriod, initial_pos_km, intial_vel_kms, timeStep):
        super().Propagate()
        timeStep            = orbitalPeriod / 1000
        totalPropagateTime  = orbitalPeriod * 2
        currentTime         = 0

        pos_state_km_ECI    = initialPos_km_ECI
        r_norm              = np.linalg.norm(pos_state_km_ECI)
        vel_state_kms_ECI   = initialVel_kms_ECI
        acc_state_kms2_ECI  = -grav_param_Earth * pos_state_km_ECI / (r_norm ** 3)

        propagatedStates.loc[len(propagatedStates)] = [currentTime, *pos_state_km_ECI.T, *vel_state_kms_ECI.T, *acc_state_kms2_ECI.T]

        while currentTime < totalPropagateTime:
            r_norm = np.linalg.norm(pos_state_km_ECI)
            # Determine new states
            [pos_state_km_ECI, vel_state_kms_ECI, acc_state_kms2_ECI] = self.EulerSolver(pos_state_km_ECI, vel_state_kms_ECI, acc_state_kms2_ECI, timeStep, grav_param_Earth)
            currentTime         += timeStep
            # Store new states in dataframe
            self.propagatedStates.loc[len(self.propagatedStates)] = [currentTime, *pos_state_km_ECI.T, *vel_state_kms_ECI.T, *acc_state_kms2_ECI.T]
            print(f"Propagating at t = {currentTime:.2f} seconds")

        return propagatedStates