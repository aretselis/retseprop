import numpy as np


def shadow_checker(x_sc, y_sc, z_sc, x_sun, y_sun, z_sun, planet_radius):
    # Computes if a spacecraft is inside a planet's shadow
    # Based on the Traditional Shadow Analysis presented in Vallado (2013)
    # Input:
    # Spacecraft position vector (x_sc, y_sc, z_sc) [meters]
    # Sun position vector (x_sun, y_sun, z_sun) [meters]
    # Output:
    # 0, if the spacecraft is in planet's shadow (umbra)
    # 1, if the spacecraft is in sunlight

    # Compute distance vector magnitudes
    r_spacecraft = np.sqrt(pow(x_sc, 2) + pow(y_sc, 2) + pow(z_sc, 2))
    r_sun = np.sqrt(pow(x_sun, 2) + pow(y_sun, 2) + pow(z_sun, 2))
    # Compute theta values
    theta = np.arccos((x_sc*x_sun+y_sc*y_sun+z_sc*z_sun)/(r_spacecraft*r_sun))
    theta_1 = np.arccos(planet_radius/r_spacecraft)
    theta_2 = np.arccos(planet_radius/r_sun)
    # Decide if we are inside or outside Earth's shadow
    if theta_1+theta_2 <= theta:
        result = 0
    else:
        result = 1
    return result
