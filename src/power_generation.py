import numpy as np
import matplotlib.pyplot as plt
import propagators
import utilities
from tqdm import trange

# Constants definition
G = 6.6743 * pow(10, -11)
M_earth = 5.977 * pow(10, 24)
mu_earth = G * M_earth
R_earth = 6378140  # [meters]

# Orbit definition
H = 500000  # Altitude, [meters]
a = R_earth + H  # Semi-major axis, [meters]
i = 97.4065  # Inclination, [deg]
e = 3.73173 * pow(10, -16)  # Eccentricity, []
Omega = 345  # Longitude of the Ascending Node [deg]
omega = 0  # Argument of pericenter [deg]
M = 0  # Mean anomaly, [deg]
T = 2 * np.pi * np.sqrt(pow(a, 3) / mu_earth)  # Orbital period, [seconds]

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
end_time = 604800
time_step = 60
x_sc_vector, y_sc_vector, z_sc_vector, vx_sc_vector, vy_sc_vector, vz_sc_vector, t_vector = \
    propagators.runge_kutta_4(x, y, z, vx, vy, vz, mu_earth, start_time, end_time, time_step)

# Compute sun position
x_sun_vector = [149597870700.0]
y_sun_vector = [0.0]
z_sun_vector = [0.0]
counter = 0
for time in trange(start_time, end_time, time_step, desc="Sun position"):
    x_sun, y_sun = utilities.sun_position_calculator(x_sun_vector[counter], y_sun_vector[counter], time)
    x_sun_vector.append(x_sun)
    y_sun_vector.append(y_sun)
    z_sun_vector.append(0.0)
    counter += 1

# Solar panel area definition
solar_panel_area = 0.0210
solar_panel_efficiency = 0.275

# Compute generated power
solar_flux = 1358
power = [0] * np.size(t_vector)
for i in trange(0, np.size(t_vector), desc="Power"):
    sunlight = utilities.shadow_checker(x_sc_vector[i], y_sc_vector[i], z_sc_vector[i], x_sun_vector[i],
                                        y_sun_vector[i], z_sun_vector[i], R_earth)
    power[i] = sunlight*solar_flux*solar_panel_area*solar_panel_efficiency

# Plot results
plt.figure()
plt.plot(t_vector, power)
plt.xlabel("Time, t, [sec]")
plt.ylabel("Power Generated, $P_{gen}$, [W]")
plt.grid()
plt.show()
