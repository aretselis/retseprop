from tqdm import tqdm

def runge_kutta_4(x_0, y_0, z_0, vx_0, vy_0, vz_0, mu, t_start, t_end, h):
    # 4th order Runge Kutta orbit propagator with no perturbations
    # Input:
    # x_0, y_0, z_0, namely the coordinates of the spacecraft's position vector [meters]
    # vx_0, vy_0, vz_0, namely the coordinates of the velocity vector [m/s]
    # mu, mu constant
    # t_start, initial time [seconds]
    # t_end, propagation end time [seconds]
    # h, integration step size [seconds]
    # Output is 7 vectors containing position, velocity and time information at each integration step

    # Log initial values
    tn = [t_start]
    xn = [x_0]
    yn = [y_0]
    zn = [z_0]
    vxn = [vx_0]
    vyn = [vy_0]
    vzn = [vz_0]
    counter = 0
    # Main RK4 loop
    progress_bar = tqdm(total=t_end, desc="Propagation")
    while tn[counter] < t_end:
        # Calculate k1 values
        k1_x = fx(vx_0)
        k1_y = fy(vy_0)
        k1_z = fz(vz_0)
        k1_vx = fv_x(x_0, y_0, z_0, mu)
        k1_vy = fv_y(x_0, y_0, z_0, mu)
        k1_vz = fv_z(x_0, y_0, z_0, mu)
        # Calculate midpoint values
        mid_x = x_0 + k1_x * h/2
        mid_y = y_0 + k1_y * h/2
        mid_z = z_0 + k1_z * h/2
        mid_vx = vx_0 + k1_vx * h/2
        mid_vy = vy_0 + k1_vy * h/2
        mid_vz = vz_0 + k1_vz * h/2
        # Calculate k2 values
        k2_x = fx(mid_vx)
        k2_y = fy(mid_vy)
        k2_z = fz(mid_vz)
        k2_vx = fv_x(mid_x, mid_y, mid_z, mu)
        k2_vy = fv_y(mid_x, mid_y, mid_z, mu)
        k2_vz = fv_z(mid_x, mid_y, mid_z, mu)
        # Calculate next midpoint values
        mid_x = x_0 + k2_x * h / 2
        mid_y = y_0 + k2_y * h / 2
        mid_z = z_0 + k2_z * h / 2
        mid_vx = vx_0 + k2_vx * h / 2
        mid_vy = vy_0 + k2_vy * h / 2
        mid_vz = vz_0 + k2_vz * h / 2
        # Calculate k3 values
        k3_x = fx(mid_vx)
        k3_y = fy(mid_vy)
        k3_z = fz(mid_vz)
        k3_vx = fv_x(mid_x, mid_y, mid_z, mu)
        k3_vy = fv_y(mid_x, mid_y, mid_z, mu)
        k3_vz = fv_z(mid_x, mid_y, mid_z, mu)
        # Calculate next midpoint values
        mid_x = x_0 + k3_x * h
        mid_y = y_0 + k3_y * h
        mid_z = z_0 + k3_z * h
        mid_vx = vx_0 + k3_vx * h
        mid_vy = vy_0 + k3_vy * h
        mid_vz = vz_0 + k3_vz * h
        # Calculate k4 values
        k4_x = fx(mid_vx)
        k4_y = fy(mid_vy)
        k4_z = fz(mid_vz)
        k4_vx = fv_x(mid_x, mid_y, mid_z, mu)
        k4_vy = fv_y(mid_x, mid_y, mid_z, mu)
        k4_vz = fv_z(mid_x, mid_y, mid_z, mu)
        # Compute r, v values and append to list
        xn.append(xn[counter] + (h / 6) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x))
        yn.append(yn[counter] + (h / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y))
        zn.append(zn[counter] + (h / 6) * (k1_z + 2 * k2_z + 2 * k3_z + k4_z))
        vxn.append(vxn[counter] + (h / 6) * (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx))
        vyn.append(vyn[counter] + (h / 6) * (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy))
        vzn.append(vzn[counter] + (h / 6) * (k1_vz + 2 * k2_vz + 2 * k3_vz + k4_vz))
        tn.append(tn[counter] + h)
        # Prepare for the next iteration and reset values
        counter += 1
        x_0 = xn[counter]
        y_0 = yn[counter]
        z_0 = zn[counter]
        vx_0 = vxn[counter]
        vy_0 = vyn[counter]
        vz_0 = vzn[counter]
        progress_bar.update(tn[counter] - progress_bar.n)
    progress_bar.close()
    return xn, yn, zn, vxn, vyn, vzn, tn


def fx(v_x):
    # Assuming x'(t)=v_x
    return v_x


def fy(v_y):
    # Assuming y'(t)=v_y
    return v_y


def fz(v_z):
    # Assuming z'(t)=v_z
    return v_z


def fv_x(x, y, z, mu):
    # Assuming v_x'(t)=-mu*x/(sqrt(x^2+y^2+z^2))^3
    return -mu*x/pow(pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1/2), 3)


def fv_y(x, y, z, mu):
    # Assuming v_y'(t)=-mu*y/(sqrt(x^2+y^2+z^2))^3
    return -mu*y/pow(pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1/2), 3)


def fv_z(x, y, z, mu):
    # Assuming v_z'(t)=-mu*z/(sqrt(x^2+y^2+z^2))^3
    return -mu*z/pow(pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1/2), 3)

