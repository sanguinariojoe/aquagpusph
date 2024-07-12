#########################################################################
#                                                                       #
#            #    ##   #  #   #                           #             #
#           # #  #  #  #  #  # #                          #             #
#          ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###           #
#          #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #          #
#          #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #          #
#          #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #          #
#                                    # #             #                  #
#                                  ##  #             #                  #
#                                                                       #
#########################################################################
#
#  This file is part of AQUA-gpusph, a free CFD program based on SPH.
#  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
#
#  AQUA-gpusph is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  AQUA-gpusph is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
#
#########################################################################

import numpy as np
import os.path as path
import aquagpusph as aqua


# Mechanical model
D = np.float32(0.062)  # Tank width
g = np.float32(9.81)   # Gravity acceleration [m/s2]
I0 = np.float32(26.9)   # Polar moment of inertia [kg m2]
m = np.float32(4.978)  # Moving mass weight [kg]
Sg = np.float32(-29.2)  # Satic moment of the rigid system [kg m]
Kdf = np.float32(0.540)  # Dry friction damping coefficient 
Bphi = np.float32(0.326)  # Linear damping coefficient

# Experimental data
FNAME = path.join('@EXAMPLE_DEST_DIR@', 'T_1-94_A100mm_water.dat')
T,XI,DXI,DDXI,THETA,DTHETA,DDTHETA = np.loadtxt(
    FNAME, delimiter=' ', skiprows=2, usecols=list(range(1, 8)), unpack=True)

# Starting Tank motion data
Theta = np.float32(0.0)
dTheta = np.float32(0.0)
ddTheta = np.float32(0.0)
DT = np.float32(0.0)

# Moment filtering data
ddtheta_list = []
dt_list = []

# Open the output file
F = open('Motion.dat', 'w')


def damping(dTheta):
    """ Compute the tank structural damping moment

    Params:
        dTheta (float): Angle time derivative dTheta/dt

    Returns:
        m (float): Tank structural damping moment
    """
    return -Kdf * np.sign(dTheta) - Bphi * dTheta


def angularForce(M, xi, dXi, theta, dTheta):
    """ Compute the tank angular acceleration

    Params:
        xi (float): Mass position
        dXi (float): Mass position time derivative dXi/dt
        theta (float): Tank angle
        dTheta (float): Tank angle time derivative dTheta/dt
    Returns:
        ddTheta (float): Tank angular acceleration ddTheta/ddt.
    """
    Mdamp = damping(dTheta)
    k0 = I0 + m*xi*xi  # ddTheta term
    k1 = 2.0*m*xi*dXi  # dTheta term
    k2 = -g*Sg         # sin(theta) term
    k3 = m*g*xi        # Independent term
    return (Mdamp + M - k1 * dTheta - k2 * np.sin(theta) - k3 * np.cos(theta)) / k0


def predictor(dx, y, dy, ddy):
    """ Performs the Heun's integrator predictor

    Params:
        dx (float): Integration step
        y (float): Current function value
        dy (float): Current function derivative dy/dx
        ddy (float): Current function second derivative ddy/ddx

    Returns:
        y (float) : Integrated function value
        dy (float): Integrated function derivative dy/dx
    """
    dy += dx * ddy
    y += dx * dy + 0.5 * dx * dx * ddy
    return y, dy


def corrector(dx, dy, ddy, ddy_in):
    """ Performs the Heun's integrator corrector

    Params:
        dx (float): Integration step
        y (float): Current function value
        dy (float): Current function derivative dy/dx
        ddy (float): Current function second derivative ddy/ddx
        ddy_in (float): Second derivative from the previous step

    Returns:
        dy (float): Integrated function derivative dy/dx
    """
    dy += 0.5 * DT * (ddy - ddy_in);
    return dy


def win_filter(m, dt, T=0.01):
    """Low-pass filter

    Params:
        m (list): Collected fluid moments
        dt (list): Collected time steps
        T (float): Filtering window

    Returns:
        m (list): Remaining moments (after droping the data out of window)
        dt (list): Remaining time steps (after droping the data out of window)
        M (float): Filtered moment value
    """
    # Drop the data out of window
    while np.sum(dt) > T:
        m = m[:-1]
        dt = dt[:-1]
    # Carry out the convolution
    f = np.exp(-(2.0 * np.cumsum(dt) / T)**2)
    M = np.sum(f * dt * m) / np.sum(f * dt)
    return m, dt, M


def main():
    global DT, Theta, dTheta, ddTheta, F, dt_list, ddtheta_list
    ddTheta_in = ddTheta  # Backup for the predictor-corrector
    # Predictor
    Theta, dTheta = predictor(DT, Theta, dTheta, ddTheta)
    # Get the new time and time step
    t = aqua.get("t")
    dt = aqua.get("dt")
    # Get the new mass position, and the expected angle
    xi = np.interp(t, T, XI)
    dXi = np.interp(t, T, DXI)
    exp_theta = np.interp(t, T, THETA)
    # Compute the angular acceleration, ddTheta
    M = aqua.get("forces_M")[0]
    dt_list.insert(0, dt)
    ddtheta_list.insert(0, angularForce(M, xi, dXi, Theta, dTheta))
    ddtheta_list, dt_list, ddTheta = win_filter(ddtheta_list, dt_list)
    # Corrector
    dTheta = corrector(dt, dTheta, ddTheta, ddTheta_in)
    DT = dt  # For the next predictor
    # Send the data to SPH
    a = np.zeros(4, dtype=np.float32)
    a[0] = Theta
    aqua.set("motion_a", a)
    dadt = np.zeros(4, dtype=np.float32)
    dadt[0] = dTheta
    aqua.set("motion_dadt", dadt)
    # Write output
    F.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        t, xi, exp_theta,
        np.degrees(a[0]), np.degrees(dadt[0]),
        M))
    F.flush()
    return True
