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


def eccentric_anomaly_calculator(time, n, tau, e):
    # Compute E using Newton-Raphson
    # Input:
    # time, [sec]
    # n, mean motion [rad/sec]
    # tau, time of pericenter [sec]
    # e, eccentricity []
    # Output
    # E, eccentric anomaly at desired time [rad]

    # Define initial value for E_0
    M = n * (time - tau)
    if M <np.pi:
        E_0 = M - e
    else:
        E_0 = M + e
    # Define f and f dot
    f = lambda E: M - E + e*np.sin(E)
    fdot = lambda E: -1 + e*np.cos(E)
    # Stopping criteria
    N = 15  # Number of significant digits to be computed
    max_repetitions = 1000000
    es = 0.5 * pow(10, (2 - N))  # Scarborough Criterion
    ea = 100
    E_prev = E_0
    repetitions = 0
    # Main Newton-Raphson loop
    while ea > es:
        repetitions = repetitions + 1
        E_next = E_prev - (f(E_prev) / fdot(E_prev))
        # if E_next == 0:
        #     return E_next
        ea = np.fabs((E_next - E_prev) * 100 / E_next)
        E_prev = E_next
        if repetitions > max_repetitions:
            raise StopIteration("Max repetitions reached without achieving desired accuracy for E!")
    E = E_next
    return E
