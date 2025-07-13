# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 15:17:09 2025

@author: waghr
"""
from Simulations.base import Simulation
import numpy as np

class PreliminaryOD(Simulation):
    def __init__(self, sat_name, sim_type):
        super().__init__('none','PreliminaryOD')
        self.sim_type = sim_type
    def Propagate(self):
        orbitalElements, orbitalPeriod = self.run()
        # timeStep is [100] steps per orbital period
        timeStep = orbitalPeriod/100
        # totalTimeofFlight is the total amount of time propagated [2] orbitalPeriod(s)
        totalTimeOfFlight = orbitalPeriod * 2
        propagatedStates = []
        for t in np.arange(0, totalTimeOfFlight + timeStep, timeStep):
            # Determine previous states, determine perturbations, propagate

            print(f"Propagating at t = {t:.2f} seconds")
        return propagatedStates