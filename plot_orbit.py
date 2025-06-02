import numpy as np
import matplotlib.pyplot as plt
"""
Inputs: 
    e: eccentricity of orbit (e = 0 circular, 0 < e < 1 elliptical, e = 1 parabolic, e > 1 hyperbolic)
    a: semi-major axis of orbit
Output:
    Plots orbit with

"""
def plot_orbit(e,a,r_vec,v_vec):
    # Determine focus 
    rp = a*(1-e)

    # True anomaly spread and finding true anomaly of r0 
    f = np.linspace(0,2 * np.pi, 1000)
    r = a*(1-e**2)/(1 + e*np.cos(f))
    #f0 = np.arccos((a * (1 - e**2) - r0) / (r0 * e))

    # Convert from polar to cartesian
    x = r*np.cos(f)
    y = r*np.sin(f)
    #x0 = r0*np.cos(f0)
    #y0 = r0*np.sin(f0)

    # Plot orbit
    plt.figure(figsize=(6,6))
    plt.plot(x, y, label="Orbit")
    plt.plot([-e*a], [0], 'ro', label="Center of Ellipse")

    # Plotting special points
    for i in range(r_vec.size):
        plt.plot()
    #plt.plot([x0],[y0], 'bo', label="Current position")
    #plt.plot([x0,0], [y0,0], linestyle="--", label="Focus to current position")
    plt.plot([0], [0], 'ko', label="Focus of Ellipse (i.e. Earth)")
    plt.plot()

    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.title(f"Orbit (e = {e}, rp = {rp})")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.show()