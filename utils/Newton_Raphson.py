import numpy as np
import matplotlib.pyplot as plt

def Newton_Raphson(M,e):
    delt = 10.0**(-10)
    Eold = 0
    Enew = M
    while abs(Enew - Eold) >= delt:
        Eold = Enew
        Enew = Enew - (Enew-e*np.sin(Enew)-M)/(1-e*np.cos(Enew))
    return(Enew)