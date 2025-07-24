import pandas as pd
import numpy as np
from utils.constants import grav_param_Earth

def EulerPropagator(orbitalPeriod, initialPos_km_ECI, initialVel_kms_ECI):
    propagatedStates = pd.DataFrame(columns=[
        'time',
        'position_I', 'position_J', 'position_K',
        'velocity_I', 'velocity_J', 'velocity_K',
        'semi_major_axis', 'eccentricity', 'true_anomaly', 'inclination', 'raan', 'arg_periapsis'
    ])
    timeStep = orbitalPeriod / 1000
    totalPropagateTime = orbitalPeriod * 2
    currentTime = 0

    pos_state_km_ECI = initialPos_km_ECI
    r_norm = np.linalg.norm(pos_state_km_ECI)
    vel_state_kms_ECI = initialVel_kms_ECI
    acc_state_kms2_ECI = -grav_param_Earth * pos_state_km_ECI / (r_norm ** 3)

    propagatedStates.loc[len(propagatedStates)] = [currentTime, *pos_state_km_ECI.T, *vel_state_kms_ECI.T, *acc_state_kms2_ECI.T]
    currentTime += timeStep

    while currentTime < totalPropagateTime:
        r_norm = np.linalg.norm(pos_state_km_ECI)
        # Determine new states
        acc_state_kms2_ECI = -grav_param_Earth * pos_state_km_ECI / (r_norm ** 3)
        vel_state_kms_ECI = vel_state_kms_ECI + acc_state_kms2_ECI * timeStep
        pos_state_km_ECI = pos_state_km_ECI + vel_state_kms_ECI * timeStep
        currentTime += timeStep
        # Store new states in dataframe
        propagatedStates.loc[len(propagatedStates)] = [currentTime, *pos_state_km_ECI.T, *vel_state_kms_ECI.T,
                                             *acc_state_kms2_ECI.T]
        print(f"Propagating at t = {currentTime:.2f} seconds")

    return propagatedStates