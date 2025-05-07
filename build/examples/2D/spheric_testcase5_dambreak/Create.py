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
hfac = 3.0
cs = 45.0
courant = 0.1
refd = 998.0
alpha = 0.0
delta = 10.0
visc_dyn = 0.000894
# Tank dimensions
h0 = 15e-2  # Reservoir water level
h1 = 38e-3  # Wetted bed level
l0 = 38e-2  # Reservoir length
l1 = 955e-2  # Wetted bed length
H = 2.0 * h0  # Total tank height
# Stimated required number of fluid particles (taking into account just the
# reservoir)
n = 50000

Vol = l0 * h0
dv = Vol / n
dr = dv**(1.0 / 2.0)

visc_dyn = max(alpha / 8.0 * refd * hfac * dr * cs, visc_dyn)

def writeParticle(output, p, n=(0.0,0.0), u=(0.0,0.0),
                  dudt=(0.0,0.0), rho=(refd), drhodt=0.0,
                  m=refd * dr**2, imove=1):
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        p[0], p[1],
        n[0], n[1],
        u[0], u[1],
        dudt[0], dudt[1],
        rho,
        drhodt,
        m,
        imove)
    output.write(string)

print("Opening output file...")
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

N = 0
string = """
    Writing reservoir fluid particles...
"""
print(string)
x = -0.5 * dr
xmin = x
while x > -l0:
    y = 0.5 * dr
    while y < h0:
        press = refd * g * (h0 - y)
        rho = refd + press / cs**2
        writeParticle(output, (x, y), rho=rho)
        N += 1
        y += dr
    x -= dr
    xmin = x
xmin += dr
print('{} particles'.format(N))

N_back = N
string = """
    Writing wetted bed fluid particles...
"""
print(string)
x = 0.5 * dr
xmax = x
while x < l1:
    y = 0.5 * dr
    while y < h1:
        press = refd * g * (h0 - y)
        rho = refd + press / cs**2
        writeParticle(output, (x, y), rho=rho)
        N += 1
        y += dr
    x += dr
    xmax = x
xmax -= dr
print('{} particles'.format(N - N_back))

N_back = N
string = """
    Writing boundary elements...
"""
print(string)
x = xmin
y = 0.0
while x < l1:
    press = refd * g * (h0 - y)
    rho = refd + press / cs**2
    writeParticle(output, (x, y), n=(0.0, -1.0), rho=rho, m=dr, imove=-3)
    N += 1
    x += dr
x = xmin - 0.5 * dr
y = 0.5 * dr
while y < H:
    press = max(0.0, refd * g * (h0 - y))
    rho = refd + press / cs**2
    writeParticle(output, (x, y), n=(-1.0, 0.0), rho=rho, m=dr, imove=-3)
    N += 1
    y += dr
x = xmax + 0.5 * dr
y = 0.5 * dr
while y < H:
    press = max(0.0, refd * g * (h0 - y))
    rho = refd + press / cs**2
    writeParticle(output, (x, y), n=(1.0, 0.0), rho=rho, m=dr, imove=-3)
    N += 1
    y += dr
print('{} particles'.format(N - N_back))

domain_min = (-1.5 * l0, -0.5 * H)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (1.5 * l1, 2.0 * H)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'N_SENSORS':str(1), 'N':str(N)}
utils.configure(data, os.path.join(script_folder, "templates"))
