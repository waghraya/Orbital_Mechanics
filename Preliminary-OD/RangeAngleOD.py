import numpy as np
from Newton_Raphson import Newton_Raphson as NR
from SEZ2ECI import SEZ2ECI as S2E
from PQW2ECI import PQW2ECI as P2E
from PQW2ECI import ECI2PQW as E2P
from plot_orbit import plot_orbit as plt_o
"""
Inputs:
    phi = 
    lam = 
    rho = scalar distance from ground station to satellite
    z = altitude of ground station relative to sea level (km)
    sig = 
    beta = 
    rho_dot = rate of change of rho with respect to time (km/s)
    sig_dot = rate change of sig with respect to time (deg/s)
    beta_dot = rate change of beta with respect to time (deg/s)
"""
    # ---------inputs-----------
phi = np.deg2rad(32.248814)
lam = np.deg2rad(-74.99)
rho = 822 #km
sig = np.deg2rad(18.0912)
beta = np.deg2rad(61.7066)
rho_dot = 3.48499169 #km/s
sig_dot = np.deg2rad(0.269604966)
beta_dot = np.deg2rad(-0.4321605433)
tof = (32)*60 #seconds

    # ---------constants----------
omega_E = np.array([[0],[0],[7.2927*10**(-5)]]) #rad/s in ECI frame
J2 = 1.08264e-3 #perturbation
rE = 6378 #km
mu = 398600 #km^3/s^2
    #---------------------[finding r1,v1]---------------------------
rho_SEZ = rho*np.array([
                        [-np.cos(sig)*np.cos(beta)],
                        [np.cos(sig)*np.sin(beta)],
                        [np.sin(sig)]
                        ])
rho_ECI = S2E(rho_SEZ,lam,phi)
rho_dot_SEZ_SEZ = rho_dot*np.array([
                        [-np.cos(sig)*np.cos(beta)],
                        [np.cos(sig)*np.sin(beta)],
                        [np.sin(sig)]
                         ]) + rho*np.array([
                             [(sig_dot*(np.sin(sig)*np.cos(beta)))+(beta_dot*np.cos(sig)*np.sin(beta))],
                             [(-sig_dot*np.sin(sig)*(np.sin(beta)))+(beta_dot*np.cos(sig)*np.cos(beta))],
                             [(sig_dot*np.cos(sig))]
                         ])
rho_dot_SEZ_ECI = S2E(rho_dot_SEZ_SEZ,lam,phi)
r_site_SEZ = np.array([[0],[0],[rE]])
r_site_ECI = S2E(r_site_SEZ,lam,phi)
r1 = r_site_ECI + rho_ECI
v1 = rho_dot_SEZ_ECI + np.cross(omega_E.T,r1.T).T

    #---------------------[r1,v1 -> OE1]----------------------------
K = np.array([[0],[0],[1]])
h = np.cross(r1.T,v1.T).T
h_hat = h/np.linalg.norm(h)
r1_hat = r1/np.linalg.norm(r1)
eps = ((np.linalg.norm(v1)**2)/2) - mu/(np.linalg.norm(r1))
n = np.cross(K.T,h_hat.T).T/(np.linalg.norm(np.cross(K.T,h_hat.T).T))
e = (np.cross(v1.T,h.T).T/mu)-r1_hat
e_hat = e/np.linalg.norm(e)
a = -mu/(2*eps)
i = float(np.acos(n.T@e_hat))

# checking domain for omega and correcting
if e_hat.T@K >= 0:
    omega = float(np.acos(n.T@e_hat))
elif e_hat.T@K < 0:
    omega = float(2*np.pi - np.acos(n.T@e_hat))
# atan2 autochecks quadrant
Ohm = float(np.atan2(n[1],n[0]))
# checking domain for f and correcting
if r1.T@v1 >= 0:
    f1 = np.acos(e_hat.T@r1_hat)
elif r1.T@v1 < 0:
    f1 = 2*np.pi - np.acos(e_hat.T@r1_hat)
E1 = 2*np.atan(np.tan(f1/2)*np.sqrt((1-np.linalg.norm(e))/(1+np.linalg.norm(e))))
# Kepler's Equation
M1 = E1 - np.linalg.norm(e)*np.sin(E1)
n = np.sqrt(mu/a**3)
e = np.linalg.norm(e)
# Convert r1,v1 to PQW
r1_PQW = E2P(r1,Ohm,omega,i)
v1_PQW = E2P(v1,Ohm,omega,i)

    #--------------------[J2 Perurbation]---------------------------
# Change of RAAN
Ohm_dot = -( ((3*n*J2)/(2*( 1- e**2)**2)) * (rE/a)**2 * np.cos(i))
Ohm = Ohm + Ohm_dot*tof
# Change of argument of periapsis
omega_dot = (3*n*J2)/(4*(1-e**2)**2)*((rE/a)**2)*(5*(np.cos(i)**2) - 1)
omega = omega + omega_dot*tof
    #----------------------[OE1 -> OE2]-----------------------------

M2 = M1 + n*tof

# Newton-Raphson method used to solve for E2
E2 = NR(M2,e)
f2 = 2 * np.atan(np.tan(E2/2) * np.sqrt((1 + np.linalg.norm(e)) / (1 - np.linalg.norm(e))))
    #---------------------[OE2 -> r2,v2]----------------------------
r = a*(1-e**2)/(1+e*np.cos(f2))
# Final r and v determined in PQW frame
r2_PQW = r*np.array([[float(np.cos(f2))],[float(np.sin(f2))],[0]])
v2_PQW = np.sqrt(mu/(a*(1-e**2)))*np.array([[float(-np.sin(f2))],[e+float(np.cos(f2))],[0]])
# Final r and v converted from PQW to ECI frame
r2 = P2E(r2_PQW,Ohm,omega,i)
v2 = P2E(v2_PQW,Ohm,omega,i)
print('r1:',r1)
print('r1_pqw:',r1_PQW)
print('r2_pqw:',r2_PQW)
r_vec = np.array([r1_PQW,r2_PQW])
v_vec = np.array([v1_PQW,v2_PQW])
print('r_vec:',r_vec)
print('v_vec:',v_vec)
plt_o(e,a,r_vec,v_vec)