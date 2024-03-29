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

BUFFER_PARTICLE = "1.0 1.0,0.0 0.0,0.0 0.0,1.0,0.0,0.0,-255,\n"

# This script is just splitting the set of particles for the 2 processes
n = [0, 0]
with open("particles.dat", "r") as fin, \
     open("particles.0.dat", "w") as out0, \
     open("particles.1.dat", "w") as out1:
    outs = [out0, out1]
    lines = fin.readlines()
    # The header
    for line in lines[:18]:
        out0.write(line)
        out1.write(line)
    # The good particles
    for line in lines[18:]:
        if line.strip() == "":
            continue
        x = float(line.split(' ')[0])
        i = 0 if x <= 0.0 else 1
        outs[i].write(line)
        n[i] += 1
    # The buffer particles
    for i in range(n[0]):
        out1.write(BUFFER_PARTICLE)
    for i in range(n[1]):
        out0.write(BUFFER_PARTICLE)
