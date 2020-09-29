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

g = 9.81
hfac = 4.0
cs = 50.0
courant = 0.1
refd = 998.0
alpha = 0.0
delta = 2.0
visc_dyn = 0.000894
# Tank dimensions
H = 0.508
L = 0.9
# Fluid
h = 0.092
# Stimated required number of fluid particles
n = 10000

# Dimensions and number of particles readjustment
# ===============================================

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
Nx *= DeLeffeDistFactor
Ny *= DeLeffeDistFactor
N = 2 * Nx + 2 * Ny

# Particles packing boundary elements
NFSx = Nx
NFSy = DeLeffeDistFactor * ny
NFS = 2 * NFSx + 2 * NFSy

# Particles generation
# ====================
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
    mass = refd * dr**2.0
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        pos[0], pos[1],
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        refd,
        0.0,
        mass,
        imove)
    output.write(string)
print('    100%')
output.close()

output = open("Tank.dat", "w")
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

output = open("FreeSurface.dat", "w")
string = """
    Writing the free-surface...
"""
print(string)
Percentage = -1
for i in range(0, NFS):
    if Percentage != (i * 100) // NFS:
        Percentage = (i * 100) // NFS
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    # Bottom
    if(i < NFSx):
        j = i
        idx = j + 0.5
        idy = 0.0
        normal = [0.0, -1.0]
    # Roof
    elif(i < 2 * NFSx):
        j = i - NFSx
        idx = j + 0.5
        idy = NFSy
        normal = [0.0, 1.0]
    # Left
    elif(i < 2 * NFSx + NFSy):
        j = i - 2 * NFSx
        idx = 0.0
        idy = j + 0.5
        normal = [-1.0, 0.0]
    # Right
    elif(i < 2 * (NFSx + NFSy)):
        j = i - 2 * NFSx - NFSy
        idx = NFSx
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

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('particlesPacking.xml', 'Fluids.xml', 'Main.xml', 'Motion.xml',
       'Settings.xml', 'SPH.xml', 'Time.xml')

radius = math.sqrt((0.5 * L)**2 + H**2)
domain_min = (-1.1 * radius, -0.5 * L * math.sin(0.5 * math.pi))
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (1.1 * radius, 1.1 * radius)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'HFS':str(hFluid), 'N':str(n), 'NBC':str(N), 'NFS':str(NFS)}
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
