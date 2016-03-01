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
from samba.dcerpc.lsa import RefDomainList

# Input data
# ==========

g = 0.0
hfac = 4.0
cs = 50.0
courant = 0.25
refd = 1.0
L = H = 1.0
U = 1.0
Re = 1000.0
T = L / U * 50
nx = ny = 100

dr = L / nx
n = nx * ny

visc_dyn = refd * U * L / Re
p0 = 3.0 * refd * U**2

alpha = 8.0 * visc_dyn / (refd * hfac * dr * cs)
delta = 1.0 if alpha < 0.03 else 0.0

Nx = nx
Ny = ny

DeLeffeDistFactor = 1
Nx = DeLeffeDistFactor * Nx
Ny = DeLeffeDistFactor * Ny
N = 2 * Nx + 2 * Ny

# Particles generation
# ====================
header = """#############################################################
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
print(header)

string = """
    Writing fluid particles...
"""
print(string)

output = open("Fluid.dat", "w")
output.write(header)

Percentage = -1
x = -0.5 * (L - dr)
for i in range(nx):
    y = -0.5 * (H - dr)
    for j in range(ny):
        if Percentage != (i * ny + j) * 100 / n:
            Percentage = (i * ny + j) * 100 / n
            if not Percentage % 10:
                string = '    {}%'.format(Percentage)
                print(string)
        imove = 1
        press = p0
        dens = refd 
        mass = dens * dr**2.0
        string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
            x, y,
            0.0, 0.0,
            0.0, 0.0,
            0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        y += dr
    x += dr
print('    100%')
output.close()

string = """
    Writing the boundary elements...
"""
print(string)

output = open("BC.dat", "w")
output.write(header)

Percentage = -1
# Top
x = -0.5 * (L - dr / DeLeffeDistFactor)
y = 0.5 * H
for i in range(0, Nx):
    if Percentage != i * 100 / N:
        Percentage = i * 100 / N
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    press = p0
    dens = refd
    mass = dr / DeLeffeDistFactor
    imove = -3
    normal = (0.0, 1.0)
    vel = (U, 0.0)
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        x, y,
        normal[0], normal[1],
        vel[0], vel[1],
        0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
    x += dr / DeLeffeDistFactor

# Bottom
x = -0.5 * (L - dr / DeLeffeDistFactor)
y = -0.5 * H
for i in range(0, Nx):
    if Percentage != (i + Nx) * 100 / N:
        Percentage = (i + Nx) * 100 / N
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    press = p0
    dens = refd
    mass = dr / DeLeffeDistFactor
    imove = -3
    normal = (0.0, -1.0)
    vel = (0.0, 0.0)
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        x, y,
        normal[0], normal[1],
        vel[0], vel[1],
        0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
    x += dr / DeLeffeDistFactor

# Left
x = -0.5 * L
y = -0.5 * (H - dr / DeLeffeDistFactor)
for i in range(0, Ny):
    if Percentage != (i + 2 * Nx) * 100 / N:
        Percentage = (i + 2 * Nx) * 100 / N
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    press = p0
    dens = refd
    mass = dr / DeLeffeDistFactor
    imove = -3
    normal = (-1.0, 0.0)
    vel = (0.0, 0.0)
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        x, y,
        normal[0], normal[1],
        vel[0], vel[1],
        0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
    y += dr / DeLeffeDistFactor

# Right
x = 0.5 * L
y = -0.5 * (H - dr / DeLeffeDistFactor)
for i in range(0, Ny):
    if Percentage != (i + 2 * Nx + Ny) * 100 / N:
        Percentage = (i + 2 * Nx + Ny) * 100 / N
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    press = p0
    dens = refd
    mass = dr / DeLeffeDistFactor
    imove = -3
    normal = (1.0, 0.0)
    vel = (0.0, 0.0)
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        x, y,
        normal[0], normal[1],
        vel[0], vel[1],
        0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
    y += dr / DeLeffeDistFactor
print('    100%')
output.close()

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('Fluids.xml', 'BCs.xml', 'Main.xml', 'Settings.xml', 'SPH.xml', 'Time.xml')

radius = math.sqrt((0.5 * L)**2 + (0.5 * H)**2)
domain_min = (-1.1 * radius, -1.1 * radius)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (1.1 * radius, 1.1 * radius)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta),  'G':str(g), 'P0':str(p0),
        'TEND':str(T), 'NFLUID':str(n), 'NBC':str(N)}
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
