# -*- coding: utf-8 -*-
"""
Created on Sun Jun 15 23:21:34 2025

@author: waghr
"""

"""
    Entry point into running scripts
"""
import sys
import traceback
from Simulations.PreliminaryOD.GibbsOD import GibbsOD
from Simulations.PreliminaryOD.RangeAngleOD import RangeAngleOD
from Simulations.PreliminaryOD.GaussLambert import GaussLambert

def getSim(simName: str):
    if simName == "RangeAngleOD":
        return RangeAngleOD.import_config("config/RangeAngleOD_config.ini")
    elif simName == "GaussLambert":
        return GaussLambert.import_config("config/GaussLambert_config.ini")
    elif simName == "GibbsOD":
        return GibbsOD.import_config("config/GibbsOD_config.ini")
    else:
        raise ValueError(f"Unknown simulation name: {simName}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python main.py <simName>")
        sys.exit(1)

    simName = sys.argv[1]

    try:
        simulation = getSim(simName)
        simulation.run()
    except Exception as e:
        traceback.print_exc()
        print(f"[ERROR] {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()