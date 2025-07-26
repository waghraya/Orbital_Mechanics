# Orbital-Determination

----- Program Description -----

This repo currently implements 3 Preliminary OD classes (Gibbs, RangeAndAngle, and Gauss-Lambert) as a means to generate the 6 orbital elements (semi-major axis, eccentricity, inclination, true anomaly, right-ascension of the ascending nodes, and argument of periapsis). And will implement different orbit maneuvers (Hohmann transfer, bi-elliptic transfer).

Within the near future, propagating the states of the orbit will be used to generate 3d visualizations and animations of the orbit.

Above is the base of how the program function where further updates will improve the fidelity by introducing different perturbations (J2, Solar Radiation Pressure, and endo-atmospheric drag). As well as state-estimation tools like an Extended Kalman Filter (EKF).

A full roadmap for the project is under the [wiki](url) roadmap page.


----- Steps to run -----

Determine the simulation you want to run:

  -PreliminaryOD:
  
    -Ground-based range and angle measurements: **'RangeAngleOD'**
    
    -Gibbs method: **'GibbsOD'**
    
    -Gauss-Lambert method: **'GaussLambert'**
    
Edit the associated configuration **'.ini'** file under the **'config'** folder.

Use the associated strings above and run in the format below while in your local clone of the repository:

    python main.py <sim_type>

Example plotFlatOrbit output:

<img width="1795" height="884" alt="image" src="https://github.com/user-attachments/assets/8e1f85de-8984-4379-9ecb-5e1c09520f9d" />

