# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 14:20:25 2025

@author: waghr
"""

class Simulation:
    def __init__(self, sim_type):
        self.sim_type = sim_type
    def getSimType(self):
        return self.sim_type
    def run(self):
        raise NotImplementedError("Subclasses must implement this method.")