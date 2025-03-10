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


g = 9.81
# Tank dimensions
L = 1.610
D = 0.600
# Fluid reservoir
B = 0.600
H = 0.300
refd = 997.0
visc_kin = 8.9e-7
# Discretization
nx = 800
dr = B / nx
ny = int(round(H / dr))
H = ny * dr
Nx = int(round(L / dr))
L = Nx * dr
A = L - B
Ny = int(round(D / dr))
D = Ny * dr
# SPH
hfac = 4.0
h = hfac * dr
courant = 0.1
Ma = 0.1
cs = (g * H)**0.5 / Ma
# Time
T = 7.15 / (g / H)**0.5
FPS = 300


e0 = 0.0
def particle(f, r):
    global e0
    m = dr**2 * refd
    e0 += m * g * r[1]
    # p = 0
    # p = 0
    # Hydrostatic
    # p = refd * g * (H - r[1])
    # Hydrostatic with a smooth function
    x = (B + r[0]) / B
    p = refd * g * (H - r[1]) * math.cos(0.5 * math.pi * x)
    rho = refd + p / cs**2
    f.write("{}, {}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, 1\n".format(
        r[0], r[1], rho, m))

def boundary(f, r, n):
    f.write("{}, {}, {}, {}, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, -3\n".format(
        r[0], r[1], n[0], n[1], refd, dr))


def sensor(f, r):
    f.write("{}, {}, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, 0\n".format(
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
n = nx * ny + 2 * (Nx + Ny)
with open("Fluid.dat", "w") as f:
    f.write(string)
    # Fluid particles
    x = -B + 0.5 * dr
    while x < 0:
        y = 0.5 * dr
        while y < H:
            particle(f, (x, y))
            y += dr
        x += dr
    # Bottom wall
    x = -B + 0.5 * dr
    while x < A:
        boundary(f, (x, 0), (0, -1))
        x += dr
    # Top wall
    x = -B + 0.5 * dr
    while x < A:
        boundary(f, (x, D), (0, 1))
        x += dr
    # Left wall
    y = 0.5 * dr
    while y < D:
        boundary(f, (-B, y), (-1, 0))
        y += dr
    # Right wall
    y = 0.5 * dr
    while y < D:
        boundary(f, (A, y), (1, 0))
        y += dr

with open("Sensors.dat", "w") as f:
    sensor(f, (A, 0.003))
    sensor(f, (A, 0.015))
    sensor(f, (A, 0.030))
    sensor(f, (A, 0.080))

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'G':str(g), 'L':str(L), 'D':str(D), 'B':str(H), 'H':str(H),
        'REFD':str(refd), 'T':str(T), 'FPS':str(FPS), 'N':str(n),
        'VISC':str(visc_kin * refd), 'E0':str(e0), 'N_SENSORS':'4'}
utils.configure(data, os.path.join(script_folder, "templates"))
