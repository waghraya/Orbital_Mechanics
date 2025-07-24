# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 15:17:09 2025

@author: waghr
"""
from Simulations.base import Simulation

class PreliminaryOD(Simulation):
    def __init__(self, sat_name, sim_type):
        super().__init__('none','PreliminaryOD')
        self.sim_type = sim_type