using .Constants

function simple_shadow_check(x_sc, y_sc, z_sc, x_sun, y_sun, z_sun, planet_radius)
    #=
    Computes if a spacecraft is inside a planet's umbra or sunlit region
    Based on the Traditional Shadow Analysis presented in Vallado (2013)
    Input:
    Spacecraft position vector (x_sc, y_sc, z_sc) [meters]
    Sun position vector (x_sun, y_sun, z_sun) [meters]
    Output:
    0, if the spacecraft is in planet's shadow (umbra)
    1, if the spacecraft is in sunlight
    =#

    # Compute distance vector magnitudes
    r_spacecraft = sqrt(x_sc^2 + y_sc^2 + z_sc^2)
    r_sun = sqrt(x_sun^2 + y_sun^2 + z_sun^2)
    # Compute theta values
    theta = acos((x_sc*x_sun+y_sc*y_sun+z_sc*z_sun)/(r_spacecraft*r_sun))
    theta_1 = acos(planet_radius/r_spacecraft)
    theta_2 = acos(planet_radius/r_sun)
    # Decide if we are inside or outside Earth's shadow
    if theta_1+theta_2 <= theta
        result = 0
    else
        result = 1
    end
    return result
    end

function advanced_shadow_check(x_sc, y_sc, z_sc, x_sun, y_sun, z_sun, planet_radius)
    #=
    Computes if a spacecraft is in umbra, penumbra or sunlit region of an orbit
    Based on the Geometrical Shadow Analysis presented in Vallado (2013)
    Input:
    Spacecraft position vector (x_sc, y_sc, z_sc) [meters]
    Sun position vector (x_sun, y_sun, z_sun) [meters]
    Output:
    0, if the spacecraft is in planet's shadow (umbra)
    1, if the spacecraft is in sunlight
    value in (0,1), if the spacecraft is in the penumbra region
    =#

    # Compute basic umbra angles
    a_umbra = deg2rad(0.264121687)       # [deg]
    a_penumbra = deg2rad(0.269007205)    # [deg]
    # Initially assume we are in sunlight, compute magnitudes and dot product
    shadow = 1.0
    r_sc_magnitude = sqrt(^(x_sc, 2) + ^(y_sc, 2) + ^(z_sc, 2))
    r_sun_magnitude = sqrt(^(x_sun, 2) + ^(y_sun, 2) + ^(z_sun, 2))
    dot_product = x_sc*x_sun + y_sc*y_sun + z_sc*z_sun
    # Decide if we should evaluate umbra/penumbra possibility
    if dot_product < 0
        dot_product = -x_sc * x_sun - y_sc * y_sun - z_sc * z_sun
        angle_s = acos(dot_product/(r_sc_magnitude*r_sun_magnitude))
        sat_horizontal = r_sc_magnitude*cos(angle_s)
        sat_vertical = r_sc_magnitude*sin(angle_s)
        x = planet_radius/sin(a_penumbra)
        pen_vertical = tan(a_penumbra)*(x+sat_horizontal)
        if sat_vertical <= pen_vertical
            y = planet_radius/sin(a_umbra)
            umb_vertical = tan(a_umbra)*(y-sat_horizontal)
            if sat_vertical <= umb_vertical
                shadow = 0.0
                else
                # Assume that Solar Intensity decreases linearly, calculate shadow value
                slope = (1-0)/(pen_vertical-umb_vertical)
                sat_vertical = sat_vertical - umb_vertical
                shadow = slope*sat_vertical
            end
        end
    end
    return shadow
    end

function eccentric_anomaly_calculator(mean_anomaly, e)
    #=
    Compute E using Newton-Raphson
    Input:
    mean_anomaly, [rad]
    e, eccentricity []
    Output
    E, eccentric anomaly at desired time [rad]
    =#

    # Define initial value for E_0)
    if mean_anomaly < pi
        E_0 = mean_anomaly - e
    else
        E_0 = mean_anomaly + e
    end
    # Define f and f dot
    f(E) = mean_anomaly - E + e*sin(E)
    fdot(E) =  -1 + e*cos(E)
    # Stopping criteria
    N = 15  # Number of significant digits to be computed
    max_repetitions = 1000000
    es = 0.5 * 10^(2.0 - N)  # Scarborough Criterion
    ea = 100
    E_prev = E_0
    repetitions = 0
    # Main Newton-Raphson loop
    while ea > es
        repetitions = repetitions + 1
        E_next = E_prev - (f(E_prev) / fdot(E_prev))
        if E_next == 0
            return E_next
        end
        ea = abs((E_next - E_prev) * 100 / E_next)
        E_prev = E_next
        if repetitions > max_repetitions
            error("Max repetitions reached without achieving desired accuracy for E!")
        end
    end
    E = E_next
    return E
    end

function orbital_elements_to_cartesian(a, e, i, Omega, omega, M, mu)
    #=
    Computes cartesian position and velocity vector given some orbital elements
    Input:
    a [m]
    e []
    i [deg]
    Omega [deg]
    omega [deg]
    M [deg]
    mu [m^3/(kg*s^2)] (GM)
    Output:
    r_vector, v_vector
    =#

    # Convert M to radians
    M = deg2rad(M)
    # Compute E using Newton-Raphson
    E = eccentric_anomaly_calculator(M, e)
    # Compute x, xdot, y, ydot on the orbital plane
    x = a * (cos(E) - e)
    y = a * sqrt(1 - e^2) * sin(E)
    r = sqrt(x^2 + y^2)
    n = sqrt(mu / (a^3))  # Mean motion
    x_dot = -(n * a^2 / r) * sin(E)
    y_dot = (n * a^2 / r) * sqrt(1 - e^2) * cos(E)
    # Rotation Matrices definition
    Omega = deg2rad(Omega)
    omega = deg2rad(omega)
    i = deg2rad(i)
    P1 = [cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1]
    P2 = [1 0 0;0 cos(i) sin(i);0 sin(i) cos(i)]
    P3 = [cos(Omega) -sin(Omega) 0;sin(Omega) cos(Omega) 0;0 0 1]
    # Compute cartesian coordinates
    x_y_vector = [x, y, 0]
    x_y_dot_vector = [x_dot, y_dot, 0]
    r_vector = *(*(*(P3, P2), P1), x_y_vector)
    v_vector = *(*(*(P3, P2), P1), x_y_dot_vector)
    return r_vector, v_vector
    end 

function sun_position_calculator(x_initial, y_initial, propagation_time)
    #=
    Rotates the sun around the earth in the ICRF frame
    Input is the initial position vector of the sun [meters] and propagation time [seconds]
    Output is the position vector of the sun [meters] after propagation
    =#
    n_dot = 1  # [deg/day]
    propagation_time = propagation_time / (60 * 60 * 24)
    theta = deg2rad(n_dot*propagation_time)
    r = sqrt(x_initial^2 + y_initial^2)
    x_final = r*cos(theta)
    y_final = r*sin(theta)
    return x_final, y_final
    end


function inclination_for_circular_sso(altitude, planet_radius=EARTH_RADIUS, planet_mu=EARTH_GRAV_CONST, planet_J2=EARTH_J2, planet_Omega_sun_syn=EARTH_OMEGA_SUN_SYN)
    #=
    Computes the inclination for a circular Sun synchronous orbit at a specific altitude
    Input: 
    altitude, (Float64), Altitude above the surface of the earth in [km]
    planet_X values default to the Earth parameters if not defined.
    planet_radius, (Float64), Radius of the planet in [m]
    planet_mu, (Float64), Gravitational parameter of the planet (G*M) in [m^3 s^−2]
    planet_J2, (Float64), 2nd order zonal harmonics coefficient in the gravitational potential expansion in []
    planet_Omega_sun_syn, (Float64), secular perturbation for RAAN as a result of J2 perurbation for a year in [rad/sec]
    Output:
    inclination, (Float64), Inclination value of the desired SSO orbit in [deg]
    =#
    altitude = altitude*1000
    semi_major_axis = planet_radius + altitude
    inclination = acosd((-2*semi_major_axis^(7/2)*planet_Omega_sun_syn)/(3*planet_radius^(2)*planet_J2*sqrt(planet_mu)))
    return inclination
    end
