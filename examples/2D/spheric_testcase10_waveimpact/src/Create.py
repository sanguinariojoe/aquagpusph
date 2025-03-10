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
hfac = 4.0
cs = 45.0
courant = 0.1
refd = 998.0
alpha = 0.0
delta = 0.1
visc_dyn = 0.000894
# Tank dimensions
H = 0.508
L = 0.9
# Fluid
h = 0.093
# Stimated required number of fluid particles
n = 10000

Vol = L * h
dv = Vol / n
dr = dv**(1.0 / 2.0)

nx = int(round(L / dr))
ny = int(round(h / dr))
n = nx * ny

hFluid = ny * dr
visc_dyn = max(alpha / 8.0 * refd * hfac * dr * cs, visc_dyn)

# Solid boundary elements
Nx = nx
Ny = int(round(H / dr)) + 1

# Correct tank dimensions
L = Nx * dr
H = Ny * dr

DeLeffeDistFactor = 1
Nx = DeLeffeDistFactor * Nx
Ny = DeLeffeDistFactor * Ny
N = 2 * Nx + 2 * Ny

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

string = """
    Writing fluid particles...
"""
print(string)

Percentage = -1
for i in range(0, n):
    if Percentage != (i * 100) // n:
        Percentage = (i * 100) // n
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    j = i
    idx = j % nx
    idy = j // nx
    imove = 1
    pos = (idx * dr - 0.5 * (L - dr),
           idy * dr + 0.5 * dr)
    press = refd * g * (hFluid - pos[1])
    dens = refd + press / cs**2 
    mass = dens * dr**2.0
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        pos[0], pos[1],
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
print('    100%')

string = """
    Writing the boundary elements...
"""
print(string)
Percentage = -1
for i in range(0, N):
    if Percentage != (i * 100) // N:
        Percentage = (i * 100) // N
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    # Bottom
    if(i < Nx):
        j = i
        idx = j + 0.5
        idy = 0.0
        normal = [0.0, -1.0]
    # Roof
    elif(i < 2 * Nx):
        j = i - Nx
        idx = j + 0.5
        idy = Ny
        normal = [0.0, 1.0]
    # Left
    elif(i < 2 * Nx + Ny):
        j = i - 2 * Nx
        idx = 0.0
        idy = j + 0.5
        normal = [-1.0, 0.0]
    # Right
    elif(i < 2 * (Nx + Ny)):
        j = i - 2 * Nx - Ny
        idx = Nx
        idy = j + 0.5
        normal = [1.0, 0.0]
    pos = (idx * dr / DeLeffeDistFactor - 0.5 * L,
           idy * dr / DeLeffeDistFactor)
    if pos[1] <= hFluid:
        press = refd * g * (hFluid - pos[1])
        dens = refd + press / cs**2 
    else:
        dens = refd
        press = 0.0
    imove = -3
    mass = dr / DeLeffeDistFactor
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        pos[0], pos[1],
        normal[0], normal[1],
        0.0, 0.0,
        0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
print('    100%')
output.close()

output = open("Sensors.dat", "w")
pos = (-0.5 * L, 0.093)
normal = (-1.0, 0.0)
dens = refd
mass = 0.0
imove = 0.0
string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
    pos[0], pos[1],
    normal[0], normal[1],
    0.0, 0.0,
    0.0, 0.0,
    dens,
    0.0,
    mass,
    imove)
output.write(string)
output.close()

radius = math.sqrt((0.5 * L)**2 + H**2)
domain_min = (-1.1 * radius, -1.1 * radius)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (1.1 * radius, 1.1 * radius)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'N_SENSORS':str(1), 'N':str(n + N)}
utils.configure(data, os.path.join(script_folder, "templates"))
