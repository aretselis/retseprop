import numpy as np


def simple_shadow_check(x_sc, y_sc, z_sc, x_sun, y_sun, z_sun, planet_radius):
    # Computes if a spacecraft is inside a planet's umbra or sunlit region
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


def advanced_shadow_check(x_sc, y_sc, z_sc, x_sun, y_sun, z_sun, planet_radius):
    # Computes if a spacecraft is in umbra, penumbra or sunlit region of an orbit
    # Based on the Geometrical Shadow Analysis presented in Vallado (2013)
    # Input:
    # Spacecraft position vector (x_sc, y_sc, z_sc) [meters]
    # Sun position vector (x_sun, y_sun, z_sun) [meters]
    # Output:
    # 0, if the spacecraft is in planet's shadow (umbra)
    # 1, if the spacecraft is in sunlight
    # value in (0,1), if the spacecraft is in the penumbra region

    # Compute basic umbra angles
    a_umbra = 0.264121687       # [deg]
    a_penumbra = 0.269007205    # [deg]
    a_umbra = np.radians(a_umbra)
    a_penumbra = np.radians(a_penumbra)
    # Initially assume we are in sunlight, compute magnitudes and dot product
    shadow = 1
    r_sc_magnitude = pow(pow(x_sc, 2) + pow(y_sc, 2) + pow(z_sc, 2), 1/2)
    r_sun_magnitude = pow(pow(x_sun, 2) + pow(y_sun, 2) + pow(z_sun, 2), 1/2)
    dot_product = x_sc*x_sun + y_sc*y_sun + z_sc*z_sun
    # Decide if we should evaluate umbra/penumbra possibility
    if dot_product < 0:
        dot_product = -x_sc * x_sun - y_sc * y_sun - z_sc * z_sun
        angle_s = np.arccos(dot_product/(r_sc_magnitude*r_sun_magnitude))
        sat_horizontal = r_sc_magnitude*np.cos(angle_s)
        sat_vertical = r_sc_magnitude*np.sin(angle_s)
        x = planet_radius/np.sin(a_penumbra)
        pen_vertical = np.tan(a_penumbra)*(x+sat_horizontal)
        if sat_vertical <= pen_vertical:
            y = planet_radius/np.sin(a_umbra)
            umb_vertical = np.tan(a_umbra)*(y-sat_horizontal)
            if sat_vertical <= umb_vertical:
                shadow = 0
            else:
                # Assume that Solar Intensity decreases linearly, calculate shadow value
                slope = (1-0)/(pen_vertical-umb_vertical)
                sat_vertical = sat_vertical - umb_vertical
                shadow = slope*sat_vertical
    return shadow


def eccentric_anomaly_calculator(mean_anomaly, e):
    # Compute E using Newton-Raphson
    # Input:
    # mean_anomaly, [rad]
    # e, eccentricity []
    # Output
    # E, eccentric anomaly at desired time [rad]

    # Define initial value for E_0)
    if mean_anomaly <np.pi:
        E_0 = mean_anomaly - e
    else:
        E_0 = mean_anomaly + e
    # Define f and f dot
    f = lambda E: mean_anomaly - E + e*np.sin(E)
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
        if E_next == 0:
            return E_next
        ea = np.fabs((E_next - E_prev) * 100 / E_next)
        E_prev = E_next
        if repetitions > max_repetitions:
            raise StopIteration("Max repetitions reached without achieving desired accuracy for E!")
    E = E_next
    return E


def orbital_elements_to_cartesian(a, e, i, Omega, omega, M, mu):
    # Computes cartesian position and velocity vector given some orbital elements
    # Input:
    # a [m]
    # e []
    # i [deg]
    # Omega [deg]
    # omega [deg]
    # M [deg]
    # mu [m^3/(kg*s^2)] (GM)
    # Output:
    # r_vector, y_vector

    # Convert M to radians
    M = np.radians(M)
    # Compute E using Newton-Raphson
    E = eccentric_anomaly_calculator(M, e)
    # Compute x, xdot, y, ydot on the orbital plane
    x = a * (np.cos(E) - e)
    y = a * np.sqrt(1 - pow(e, 2)) * np.sin(E)
    r = np.sqrt(pow(x, 2) + pow(y, 2))
    n = np.sqrt(mu / pow(a, 3))  # Mean motion
    x_dot = -(n * pow(a, 2) / r) * np.sin(E)
    y_dot = (n * pow(a, 2) / r) * np.sqrt(1 - pow(e, 2)) * np.cos(E)
    # Rotation Matrices definition
    Omega = np.radians(Omega)
    omega = np.radians(omega)
    i = np.radians(i)
    P1 = np.array([[np.cos(omega), -np.sin(omega), 0],
                   [np.sin(omega), np.cos(omega), 0],
                   [0, 0, 1]])
    P2 = np.array([[1, 0, 0],
                   [0, np.cos(i), np.sin(i)],
                   [0, np.sin(i), np.cos(i)]])
    P3 = np.array([[np.cos(Omega), -np.sin(Omega), 0],
                   [np.sin(Omega), np.cos(Omega), 0],
                   [0, 0, 1]])
    # Compute cartesian coordinates
    x_y_vector = np.array([x, y, 0])
    x_y_dot_vector = np.array([x_dot, y_dot, 0])
    r_vector = np.matmul(np.matmul(np.matmul(P3, P2), P1), x_y_vector)
    v_vector = np.matmul(np.matmul(np.matmul(P3, P2), P1), x_y_dot_vector)
    return r_vector, v_vector


def sun_position_calculator(x_initial, y_initial, propagation_time):
    # Rotates the sun around the earth in the ICRF frame
    # Input is the initial position vector of the sun [meters] and propagation time [seconds]
    # Output is the position vector of the sun [meters] after propagation
    n_dot = 1  # [deg/day]
    propagation_time = propagation_time / (60 * 60 * 24)
    theta = n_dot*propagation_time
    theta = np.radians(theta)
    r = np.sqrt(pow(x_initial, 2) + pow(y_initial, 2))
    x_final = r*np.cos(theta)
    y_final = r*np.sin(theta)
    return x_final, y_final
