#! /usr/bin/env python
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
#  This file is part of AQUAgpusph, a free CFD program based on SPH.
#  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
#
#  AQUAgpusph is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  AQUAgpusph is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
#
#########################################################################

# This script is just generating the set of particles. New particles are
# not generated each time the test is executed.
import numpy as np


N = 500
r = np.zeros((N, 4), dtype=float)
r[:, 0:-1] = np.random.rand(N, 3)
# Let's create a completelly alone particle
r[0, :] -= 2
# Create a tensile instability
r[-2, :] = r[-1, :]
# Create a tensile instability outside the computational domain
r[1, :] += 2
r[2, :] = r[1, :]

with open("particles.dat", "w") as output:
    for i in range(N):
        output.write("{} {} {} 0.0\n".format(
            r[i, 0], r[i, 1], r[i, 2]))
