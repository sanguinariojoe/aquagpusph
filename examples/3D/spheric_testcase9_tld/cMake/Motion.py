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


D = np.float32(0.062)  # Tank width
g = np.float32(9.81)   # Gravity acceleration [m/s2]
I0 = np.float32(26.9)   # Polar moment of inertia [kg m2]
m = np.float32(4.978)  # Moving mass weight [kg]
Sg = np.float32(-29.2)  # Satic moment of the rigid system [kg m]
Kdf = np.float32(0.540)  # Dry friction damping coefficient 
Bphi = np.float32(0.326)  # Linear damping coefficient

FNAME = path.join('@EXAMPLE_DEST_DIR@', 'T_1-94_A100mm_water.dat')
# We must skip the column zero, which is empty due to the heading space of each
# file line
T,XI,DXI,DDXI,THETA,DTHETA,DDTHETA = np.loadtxt(
    FNAME, delimiter=' ', skiprows=2, usecols=list(range(1, 8)), unpack=True)

# Starting Tank motion data
Theta = np.float32(0.0)
dTheta = np.float32(0.0)
ddTheta = np.float32(0.0)
DT = np.float32(0.0)

# Open the output file
F = open('Motion.dat', 'w')


def damping(dTheta):
    """ Compute the tank structural damping moment.
    Params:
        dTheta: Angle time derivative dTheta/dt
    Returns:
        Tank structural damping moment.
    """
    return -Kdf * np.sign(dTheta) - Bphi * dTheta


def angularForce(M, xi, dXi, theta, dTheta):
    """ Compute the tank angular acceleration.
    Params:
        xi: Mass position
        dXi: Mass position time derivative dXi/dt
        theta: Tank angle
        dTheta: Tank angle time derivative dTheta/dt
    Returns:
        Tank angular acceleration ddTheta/dt.
    """
    Mdamp = damping(dTheta)
    k0 = I0 + m*xi*xi  # ddTheta term
    k1 = 2.0*m*xi*dXi  # dTheta term
    k2 = -g*Sg         # sin(theta) term
    k3 = m*g*xi        # Independent term
    return (Mdamp + M - k1 * dTheta - k2 * np.sin(theta) - k3 * np.cos(theta)) / k0


def predictor(dx, y, dy, ddy):
    """ Performs predictor Leap-Frog stage.
    Params:
        dx: Integration step
        y: Current function value
        dy: Current function derivative dy/dx
        ddy: Current function second derivative ddy/ddx
    """
    dy += dx * ddy
    y += dx * dy + 0.5 * dx * dx * ddy
    return y, dy


def corrector(dx, dy, ddy, ddy_in):
    """ Performs corrector Leap-Frog stage.
    Params:
        dx: Integration step
        y: Current function value
        dy: Current function derivative dy/dx
        ddy: Current function second derivative ddy/ddx
        ddy_in: Second derivative from the previous step
    """
    dy += 0.5 * DT * (ddy - ddy_in);
    return dy


def main():
    global DT, Theta, dTheta, ddTheta, F
    ddTheta_in = ddTheta  # Backup for the predictor-corrector
    # Predictor
    Theta, dTheta = predictor(DT, Theta, dTheta, ddTheta)
    # Get the new time and time step
    t = aqua.get("t")
    dt = aqua.get("dt")
    # Get new mass position (and experimental resulting angle)
    xi = np.interp(t, T, XI)
    dXi = np.interp(t, T, DXI)
    exp_theta = np.interp(t, T, THETA)
    # Compute ddTheta
    M = aqua.get("forces_M")[0]
    ddTheta = angularForce(M, xi, dXi, Theta, dTheta)
    # Corrector
    dTheta = corrector(dt, dTheta, ddTheta, ddTheta_in)
    DT = dt  # Store it for the predictor in the following time step
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