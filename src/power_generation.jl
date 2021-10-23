include("propagators.jl")
include("utilities.jl")
using Plots, ProgressBars, LaTeXStrings

function main()
    # Constants definition
    G = 6.6743 * (10^-11)
    M_earth = 5.977 * (10.0^24)
    mu_earth = G * M_earth
    R_earth = 6378140  # [meters]

    # Orbit definition
    H = 500000  # Altitude, [meters]
    a = Float64(R_earth + H)  # Semi-major axis, [meters]
    i = 97.4065  # Inclination, [deg]
    e = 3.73173 * 10^(-16)  # Eccentricity, []
    Omega = 345  # Longitude of the Ascending Node [deg]
    omega = 0  # Argument of pericenter [deg]
    M = 0  # Mean anomaly, [deg]
    T = 2 * Float64(pi) * sqrt((a^3) / mu_earth)  # Orbital period, [seconds]

    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a, e, i, Omega, omega, M,  mu_earth)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate orbit (all times in seconds)
    start_time = 0
    end_time = 3*T
    time_step = 60
    x_sc_vector, y_sc_vector, z_sc_vector, vx_sc_vector, vy_sc_vector, vz_sc_vector, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_earth, start_time, end_time, time_step)

    # Compute sun position
    x_sun_vector = [139861683376.804]
    y_sun_vector = [2359460644.848]
    z_sun_vector = [59677952894.078]
    r_vector = Vector{Float64}(undef, size(x_sc_vector, 1))
    counter = 1
    println("\nComputing Sun Position: ")
    for time in ProgressBar(start_time:time_step:end_time)
        x_sun, y_sun = sun_position_calculator(x_sun_vector[counter], y_sun_vector[counter], time)
        push!(x_sun_vector, x_sun)
        push!(y_sun_vector, y_sun)
        push!(z_sun_vector, 0.0)
        counter += 1
    end

    # Solar panel area definition
    solar_panel_area = 0.021049
    solar_panel_efficiency = 0.284

    # Compute generated power
    solar_flux = 1367
    power = Vector{Float64}(undef, size(t_vector, 1))
    println("\nComputing Power: ")
    for i in ProgressBar(1:size(t_vector, 1))
        sunlight = advanced_shadow_check(x_sc_vector[i], y_sc_vector[i], z_sc_vector[i], x_sun_vector[i], y_sun_vector[i], z_sun_vector[i], R_earth)
        power[i] = sunlight * solar_flux * solar_panel_area * solar_panel_efficiency
    end

    # Plot results
    pyplot()
    plot(t_vector, power)
    xlabel!("Time, \$t\$, [sec]")
    ylabel!("Power Generated, \$P_{gen}\$, [W]")
end

main()