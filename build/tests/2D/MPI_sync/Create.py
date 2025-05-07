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
PROCS = 2

r = np.random.rand(N, 2)

n = [0]
dn = N // PROCS
for proc in range(PROCS):
    n.append(n[-1] + dn)
n[-1] = N

for proc in range(PROCS):
    next_proc = proc + 1
    if next_proc == PROCS:
        next_proc = 0
    with open("particles_{}.dat".format(proc), "w") as output:
        for i in range(n[proc], n[proc + 1]):
            output.write("{} {} {}\n".format(
                r[i, 0], r[i, 1], next_proc))
