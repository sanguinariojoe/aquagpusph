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

# Rüchardt's experiment generalization
# See Dittman, Richard H., and Mark W. Zemansky. "Heat and thermodynamics".

import os.path as path
import math

# Input data
# ==========

# SPH resolution
ny = nx = 400
hfac = 4.0
# Box size (on still state)
L = H = 1.0
# Fluid properties
refd = 1.0
cs = 1.0
# Match number
Ma = 0.01
# Maximum expansion before we get negative pressure
rho_eq = (1. + Ma) * refd
F = cs * cs * refd * H * Ma
x_max = Ma * L
# Mass ratio between the fluid and the piston
M = (1. + Ma) * refd * L * H
mass = 0.1 * M
# Initial fluid compression
x_0 = -0.5 * x_max
V_0 = (L + x_0) * H
rho = M / V_0
# SPH vars
dr = L / nx
h = hfac * dr
courant = 0.5
# Simulation time
T0_factor = -math.log(mass / M)
T0 = T0_factor * 2. * math.pi / cs * (H * L**2 * mass / M)**0.5
T = 1.0 * T0
FPS = 100 / T0

def particle(f, r):
    m = dr**2 * rho
    f.write("{}, {}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, 1\n".format(
        r[0], r[1], rho, m))

def boundary(f, r, n, dr):
    f.write("{}, {}, {}, {}, 0.0, 0.0, 0.0, 0.0, {}, 0.0, {}, -3\n".format(
        r[0], r[1], n[0], n[1], rho, dr))

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
n = nx * ny + 2 * nx + 2 * ny
with open("Fluid.dat", "w") as f:
    f.write(string)
    dr_0 = (L + x_0) / nx
    # Fluid particles
    x = -L + 0.5 * dr_0
    while x < x_0:
        y = -0.5 * H + 0.5 * dr
        while y < 0.5 * H:
            particle(f, (x, y))
            y += dr
        x += dr_0
    # Left wall
    x = -L
    y = -0.5 * H + 0.5 * dr
    while y < 0.5 * H:
        boundary(f, (x, y), (-1.0, 0.0), dr)
        y += dr
    # Right wall
    x = x_0
    y = -0.5 * H + 0.5 * dr
    while y < 0.5 * H:
        boundary(f, (x, y), (1.0, 0.0), dr)
        y += dr
    # Bottom wall
    y = -0.5 * H
    x = -L + 0.5 * dr_0
    while x < x_0:
        boundary(f, (x, y), (0.0, -1.0), dr_0)
        x += dr_0
    # Top wall
    y = 0.5 * H
    x = -L + 0.5 * dr_0
    while x < x_0:
        boundary(f, (x, y), (0.0, 1.0), dr_0)
        x += dr_0

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('Fluids.xml', 'Main.xml', 'Settings.xml', 'SPH.xml', 'Time.xml',
       'Spring.xml', 'plot_f.py', 'plot_e.py', 'plot_x.py')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'L':str(L), 'H':str(H), 'REFD':str(refd), 'T':str(T), 'T0':str(T0),
        'FPS':str(FPS), 'NX':str(nx), 'NY':str(ny),
        'N':str(n), 'X0':str(x_0), 'F':str(F), 'M':str(mass)}
for fname in XML:
    # Read the template
    f = open(path.join(templates_path, fname), 'r')
    txt = f.read()
    f.close()
    # Replace the data
    for k in data.keys():
        txt = txt.replace('{{' + k + '}}', data[k])
    # Write the file
    f = open(fname, 'w')
    f.write(txt)
    f.close()
