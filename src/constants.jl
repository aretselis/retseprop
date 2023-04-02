module Constants

export EARTH_GRAV_CONST, EARTH_RADIUS, EARTH_J2, EARTH_J3, 
        EARTH_J4, EARTH_OMEGA_SUN_SYN

const EARTH_GRAV_CONST = 3.986004418e14 # Gravitational constant for Earth, in m^3/s^2
const EARTH_RADIUS = 6.378137e6 # Mean radius of Earth, in m
const EARTH_J2 = 1.0826267e-3 # 2nd order zonal harmonics coefficient in the gravitational potential expansion of the Earth in []
const EARTH_J3 = -0.0000025327 # 3rd order zonal harmonics coefficient in the gravitational potential expansion of the Earth in []
const EARTH_J4 = -0.0000016196 # 4th order zonal harmonics coefficient in the gravitational potential expansion of the Earth in []
const EARTH_OMEGA_SUN_SYN = 1.991063853e-7 #secular perturbation for RAAN as a result of J2 perurbation for a year for the Earth in [rad/sec]

end