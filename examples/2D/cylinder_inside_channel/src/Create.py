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

g = 0.0
hfac = 3.0
courant = 0.05
courant_ramp_iters = 1000
courant_ramp_factor = 0.001
refd = 1.0
alpha = 0.0
delta = 20.0
U = 1.0
cs = 10.0 * U
Re = 100.0
p0 = 3.0 * refd * U**2
# Cylinder and Channel dimensions
D = 1.0
L = 12.0 * D
H = 5.0 * D
# Position of the cylinder
x_cyl = 4.0 * D
y_cyl = 0.0 * D
# Number of fluid particles in y direction
ny = 300


sep = 2.0
dr = H / ny
h = hfac * dr
domain_min = (-3.0 * sep * h, -(0.5 * H + 2.0 * sep * h))
domain_max = (L + 6.0 * sep * h, 0.5 * H + 3.0 * sep * h)

n_buffer_x = max(int(8.0 * sep * hfac), int(D / dr))
n_buffer_y = ny
n_buffer = n_buffer_x * n_buffer_y

visc_dyn = refd * U * D / Re
visc_dyn = max(alpha / 8.0 * refd * hfac * dr * cs, visc_dyn)

print("Opening fluid particles output file...")
output = open("Fluid.dat", "w")
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
    Writing fluid particles...
"""
print(string)

n = 0
Percentage = -1
x = -0.5 * dr - sep * h
while x <= L + sep * h + 0.5 * dr:
    percentage = int(round((x + sep * h) / (L + 2.0 * sep * h + dr) * 100))
    if Percentage != percentage:
        Percentage = percentage
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    y = -0.5 * H + 0.5 * dr
    while y < 0.5 * H:
        # Avoid the particles inside the cylinder
        dist = math.sqrt((x - x_cyl)**2 + (y - y_cyl)**2)
        if dist < 0.5 * (D + dr):
            y += dr
            continue
        n += 1
        imove = 1
        vx = 0.0
        vy = 0.0
        press = p0
        dens = refd + (press - p0) / cs**2
        mass = refd * dr**2.0
        string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
            x, y,
            0.0, 0.0,
            vx, vy,
            0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        y += dr
    x += dr

string = """
    Writing the boundaries...
"""
print(string)

Percentage = -1
x = -0.5 * dr - sep * h
while x <= L + sep * h + 0.5 * dr:
    percentage = int(round((x + sep * h) / (L + 2.0 * sep * h + dr) * 100))
    if Percentage != percentage:
        Percentage = percentage
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    for i, y in enumerate([-0.5 * H, 0.5 * H]):
        n += 1
        imove = -3
        vx = 0.0
        vy = 0.0
        normal_x = 0
        normal_y = 1 if i else -1
        press = p0
        dens = refd + (press - p0) / cs**2 
        mass = dr
        string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
            x, y,
            normal_x, normal_y,
            vx, vy,
            0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
    x += dr

string = """
    Writing buffer particles...
"""
print(string)

Percentage = -1
x = domain_max[0] + sep * h
y = domain_max[1] + sep * h
for i in range(n_buffer):
    percentage = int(round((i + 1.0) / n_buffer * 100))
    if Percentage != percentage:
        Percentage = percentage
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    n += 1
    imove = -255
    vx = 0.0
    vy = 0.0
    press = 0.0
    dens = refd
    mass = refd * dr**2.0
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        x, y,
        0.0, 0.0,
        vx, vy,
        0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
output.close()
print('{} particles written'.format(n))

print("Opening cylinder boundary elements output file...")
output = open("Cylinder.dat", "w")
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

string = """
    Writing the boundary elements...
"""
print(string)
Percentage = -1
theta = 0.0
dtheta = dr / (0.5 * D)
n_cyl = int(round(2.0 * math.pi / dtheta))
dtheta = 2.0 * math.pi / n_cyl
n_cyl = 0  # Avoid rounding errors
while theta < 2.0 * math.pi:
    percentage = int(round(theta / (2.0 * math.pi) * 100))
    if Percentage != percentage:
        Percentage = percentage
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)

    imove = -3
    x = x_cyl + 0.5 * D * math.cos(theta)
    y = y_cyl + 0.5 * D * math.sin(theta)
    vx = 0.0
    vy = 0.0
    normal_x = -math.cos(theta)
    normal_y = -math.sin(theta)
    press = p0
    dens = refd + (press - p0) / cs**2 
    mass = dr
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        x, y,
        normal_x, normal_y,
        vx, vy,
        0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)

    n_cyl += 1
    theta += dtheta
output.close()
print('{} boundary elements written'.format(n_cyl))


domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'N':str(n), 'NY':str(ny), 'L':str(L), 'H':str(H), 'D':str(D),
        'U':str(U), 'P0':str(p0),
        'NCYL':str(n_cyl), 'XCYL':str(x_cyl), 'YCYL':str(y_cyl),
        'COURANT_RAMP_ITERS':str(courant_ramp_iters),
        'COURANT_RAMP_FACTOR':str(courant_ramp_factor)}
utils.configure(data, os.path.join(script_folder, "templates"))
