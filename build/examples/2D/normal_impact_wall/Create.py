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

import os
import sys
script_folder = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_folder, "../../"))
import aqua_example_utils as utils
import math

nx = 400
ny = nx / 2
U = 1.0
L = 1.0
H = 0.5
Ma = 0.01
refd = 1.0
dr = L / nx
hfac = 4.0
h = hfac * dr
courant = 0.5
T = 3
FPT = 50

def particle(f, r):
    uy = -U if r[1] > 0.0 else U
    m = dr**2 * refd
    f.write("{}, {}, 0.0, {}, 0.0, 0.0, {}, 0.0, {}, 1\n".format(
        r[0], r[1], uy, refd, m))

def boundary(f, r):
    f.write("{}, {}, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, -3\n".format(
        r[0], r[1], refd, dr))

string = """#############################################################
#                                                           #
#    #    ##   #  #   #                           #         #
#   # #  #  #  #  #  # #                          #         #
#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###       #
#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #      #
#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #      #
#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #      #
#                            # #             #              #
#                          ##  #             #              #
#                                                           #
#############################################################
"""
n = nx * ny + 3 * nx
with open("Fluid.dat", "w") as f:
    f.write(string)
    # Fluid particles
    x = -(nx / 2 - 0.5) * dr
    while x < nx / 2 * dr:
        y = 0.5 * dr
        while y < ny * dr:
            particle(f, (x, y))
            y += dr
        x += dr
    # Boundary particles
    x = -(3 * nx / 2 - 0.5) * dr
    while x < 3 * nx / 2 * dr:
        boundary(f, (x, 0))
        x += dr

data = {'DR':str(dr), 'HFAC':str(hfac), 'MA':str(Ma), 'COURANT':str(courant),
        'L':str(L), 'H':str(H), 'U':str(U), 'REFD':str(refd), 'T':str(T * Ma),
        'FPS':str(int(round(FPT / Ma))), 'NX':str(nx), 'NY':str(ny),
        'N':str(n), 'E0':str(0.5 * L * H * refd * U**2)}
utils.configure(data, os.path.join(script_folder, "templates"))
