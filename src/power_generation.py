import numpy as np
import propagators
import utilities


# Constants definition
G = 6.6743*pow(10, -11)
M_earth = 5.977*pow(10, 24)
mu_earth = G*M_earth
R_earth = 6378140  # [meters]

# Orbit definition
H = 500000  # Altitude, [meters]
a = R_earth + H  # Semi-major axis, [meters]
i = 97.4065  # Inclination, [deg]
e = 3.73173*pow(10, -16)  # Eccentricity, []
Omega = 345  # Longitude of the Ascending Node [deg]
omega = 0  # Argument of pericenter [deg]
M = 0  # Mean anomaly, [deg]
T = 2 * np.pi * np.sqrt(pow(a, 3)/mu_earth)  # Orbital period, [seconds]

# Compute initial position and velocity vector from orbital elements
r_vector, v_vector = utilities.orbital_elements_to_cartesian(a, e, i, Omega, omega, M, mu_earth)
x = r_vector[0]
y = r_vector[1]
z = r_vector[2]
vx = v_vector[0]
vy = v_vector[1]
vz = v_vector[2]

# Propagate orbit (all times in seconds)
start_time = 0
end_time = T
time_step = 604800
xn, yn, zn, vxn, vyn, vzn, tn = \
    propagators.runge_kutta_4(x, y, z, vx, vy, vz, mu_earth, start_time, end_time, time_step)

