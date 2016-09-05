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

# See:
# Iason Zisis, Ramzi Messahel, Abdelaziz Boudlal, Bas van der Linden and
# Barry Koren. "Validation of robust SPH schemes for fully compressible
# multiphase flows". Int. Jnl. of Multiphysics. Vol 9(3), Pp 225-234, 2015

import os.path as path
import math

# Input data
# ==========

hfac = 2.0                 # nu factor in the paper
courant = 0.01             # Useless if dt is not None
refd_1 = 1.0
refd_2 = 0.125
press_1 = 1.0
press_2 = 0.1
alpha = 0.0                # Artificial viscosity factor
delta = 0.0                # delta-SPH factor (Hashemi and Manzari model)
gamma = 5.0 / 3.0
l_1 = h_1 = 0.5
l_2 = h_2 = 2.5
L = H = 3.0
end_time = 0.5             # None to automatically compute it
dt = 1.0E-4                # None to automatically compute it
nx, ny = 100 * L, 100 * H

# Dimensions and number of particles readjustment
# ===============================================
dr = L / nx
n = nx * ny

# Resulting simulation parameters
# ===============================
cs = math.sqrt(max(gamma * press_1 / refd_1, gamma * press_2 / refd_2))
visc_dyn = alpha / 8.0 * min(refd_1, refd_2) * hfac * dr * cs
dt = dt or courant * hfac * dr / cs
end_time = end_time or 0.75 * min(L, H) / cs

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
while x <= L:
    y = 0.5 * dr
    while y <= H:
        N += 1
        imove = 1 if x < l_2 and y < h_2 else -1
        dens = refd_1 if x < l_1 and y < h_1 else refd_2
        press = press_1 if x < l_1 and y < h_1 else press_2
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
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'GAMMA':str(gamma),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'L':str(L), 'H':str(H),
        'N':str(N), 'DT':str(dt), 'END_TIME':str(end_time)}
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