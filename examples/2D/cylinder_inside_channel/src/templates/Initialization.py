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
import aquagpusph as aqua


# Outside set constants
U = {{U}}
L = {{L}} / 10.0

# We are using a Wendland 1D kernel for the acceleration profile, such that
# its integral (velocity profile) is a smooth function starting from 0 (for
# q=-2) and ending at 1 (for q=2)
Q = np.array(np.linspace(-2, 2, num=10000), dtype=np.float32)
W = 3.0 / 64.0 * (1.0 + 2.0 * abs(Q)) * (2.0 - abs(Q))**4
F = []
for i in range(0, len(Q)):
    F.append(np.trapz(W[:i+1], x=Q[:i+1]))
F = np.array(F, dtype=np.float32)


def expandTime(t):
    """Translate the time variable to the non-dimensional Q space"""
    # First we must adimensionalize the time
    T = t * U / L
    # Now we must expand it from the range (0, 1) to the range (-2, 2)
    return 4.0 * T - 2.0    


def evalDUdt(t):
    """ Compute the acceleration """
    T = expandTime(t)
    # 4 * U / L due to the variable change
    return 4.0 * U / L * np.interp(T, Q, W) * U


def evalU(t):
    """ Compute the velocity """
    T = expandTime(t)
    return np.interp(T, Q, F) * U


def main():
    # Get the time instant
    t = aqua.get("t")
    dt = aqua.get("dt")
    # Evaluate the acceleration at t + 0.5 * dt, and the velocity at t + dt
    dUdt = evalDUdt(t + 0.5 * dt)
    U = evalU(t + dt)
    # Set the modified variables
    aqua.set("U", float(U))
    aqua.set("dUdt", float(dUdt))
    aqua.set("inlet_U", float(U))
    aqua.set("outlet_U", float(U))
    return True
