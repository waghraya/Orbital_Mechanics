import numpy as np
from plot_orbit import plot_orbit as plt_o

"""
    Inputs:
        Three position vectors: r1, r2, r3        
"""

   #------inputs-----
mu = 0 #change to be inputted from .ini
r1 = np.array([[0.0],[0.0],[0.0]])
r2 = np.array([[0.0],[0.0],[0.0]])
r3 = np.array([[0.0],[0.0],[0.0]])  

r1_norm = np.linalg.norm(r1)
r2_norm = np.linalg.norm(r2)
r3_norm = np.linalg.norm(r3)

 #----Define D,N,S-----
D = np.cross(r1,r2) + np.cross(r2,r3) + np.cross(r3,r1)
N = r1_norm*np.cross(r2,r3) + r2_norm*np.cross(r3,r1) + r3_norm*np.cross(r1,r2)
S = (r2_norm-r3_norm)*r1+(r3_norm-r1_norm)*r2+(r1_norm-r2_norm)*r3

 #----Solve for v----

gibbs_intermediate = np.sqrt(mu/np.dot(N,D))
v1 = (gibbs_intermediate/r1_norm)*(np.cross(D,r1)) + gibbs_intermediate*S
v2 = (gibbs_intermediate/r2_norm)*(np.cross(D,r2)) + gibbs_intermediate*S
v3 = (gibbs_intermediate/r3_norm)*(np.cross(D,r3)) + gibbs_intermediate*S

 #----RV2OE-----
#[array of OEs] = RV2OE(r,v)