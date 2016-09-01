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

hfac = 2.0
cs = 100.0
courant = 0.1
refd_1 = 1.0
refd_2 = 0.125
press_1 = 1.0
press_2 = 0.1
alpha = 0.0
delta = 1.0
gamma = 1.4
l = h = 0.1
L = H = 0.5
end_time = 0.5
# Number of fluid particles in y direction
nx = ny = 100

# Dimensions and number of particles readjustment
# ===============================================
dr = L / nx
n = nx * ny

visc_dyn = alpha / 8.0 * min(refd_1, refd_2) * hfac * dr * cs

# Particles generation
# ====================
print("Opening output file...")
output = open("Particles.dat", "w")
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
output.write(string)
print(string)

string = """
    Writing particles...
"""
print(string)

N = 0
x = 0.5 * dr
while x < L + 2 * hfac * dr:
    y = 0.5 * dr
    while y < H + 2 * hfac * dr:
        N += 1
        imove = 1 if x < L and y < H else -1
        dens = refd_1 if x < l and y < h else refd_2
        press = press_1 if x < l and y < h else press_2
        e_int = press / ((gamma - 1.0) * dens)
        mass = dens * dr**2
        string = ("{} {}, " * 4 + "{}, {}, {}, {}, {}, {}\n").format(
            x, y,
            0.0, 0.0,
            0.0, 0.0,
            0.0, 0.0,
            dens,
            0.0,
            e_int,
            0.0,
            mass,
            1.0,        # Omega, required by the fixed particles
            imove)
        output.write(string)
        y += dr
    x += dr
print('OK')

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('BCs.xml', 'EOS.cl', 'EOS.xml', 'Fluids.xml', 'Main.xml', 'Settings.xml',
       'SPH.xml', 'Time.xml')

domain_min = (-0.5 * L, -0.5 * H)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (1.5 * L, 1.5 * H)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max,
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'L':str(L), 'H':str(H),
        'N':str(N), 'END_TIME':str(end_time), 'GAMMA':str(gamma)}
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