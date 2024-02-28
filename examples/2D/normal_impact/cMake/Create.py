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

import os.path as path
import math

# Input data
# ==========

nx = ny = 400
U = 1.0
L = H = 1.0
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
n = nx * ny
with open("Fluid.dat", "w") as f:
    f.write(string)
    x = -(nx / 2 - 0.5) * dr
    while x < nx / 2 * dr:
        y = -(ny / 2 - 0.5) * dr
        while y < ny / 2 * dr:
            particle(f, (x, y))
            y += dr
        x += dr

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('Fluids.xml', 'Main.xml', 'Settings.xml', 'SPH.xml', 'Time.xml',
       'plot_e.py')

data = {'DR':str(dr), 'HFAC':str(hfac), 'MA':str(Ma), 'COURANT':str(courant),
        'L':str(L), 'H':str(H), 'U':str(U), 'REFD':str(refd), 'T':str(T * Ma),
        'FPS':str(int(round(FPT / Ma))), 'NX':str(nx), 'NY':str(ny),
        'N':str(n), 'E0':str(0.5 * L * H * refd * U**2)}
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
