# -*- coding: utf-8 -*-
"""
Created on Sun Jun 15 23:32:58 2025

@author: waghr
"""

class OrbitDetermination:
    def __init__(self, observations):
        self.observations = observations
        
    def computeOrbit(self):
        raise NotImplementedError('Implement in subclass')