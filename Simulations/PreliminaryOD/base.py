# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 15:17:09 2025

@author: waghr
"""
from Simulations.base import Simulation

class PreliminaryOD(Simulation):
    def __init__(self, sim_type):
        self.sim_type = sim_type
    def getSimType(self):
        sim_type = 'Preliminar_OD'
        return sim_type