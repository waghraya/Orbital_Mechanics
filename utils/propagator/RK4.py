import numpy as np
from base import Propagator

class RK4(Propagator):
    def __init__(self, time_step, initial_pos_ECI, initial_vel_ECI):
        propagator_type = '4th-order Runge-Kutta'
        super().__init__(time_step, propagator_type, initial_pos_ECI, initial_vel_ECI)

    def Propagate(self, orbitalPeriod, initialPos_km_ECI, initialVel_kms_ECI):
        super().Propagate()
        #Solve K1 (Euler)
        # K1 = f(tn, yn)
        k1 = self.two_body_ODE(pos_km, vel_kms, grav_param)
        #Solve K2 (with K1)
        # K2 = f(tn + h/2, yn + h/2 * k1)
        k2 = self.two_body_ODE(pos_km, vel_kms, grav_param)

        #Solve K3 (with K2)
        # K3 = f(tn + h/2, yn + h/2 * k2)

        #Solve K4 (with K3)
        # K4 = f(tn + h, yn + h * k3)


        # yn+1 = yn + h/6*(k1 + 2*k2 + 2*k3 + k4)
