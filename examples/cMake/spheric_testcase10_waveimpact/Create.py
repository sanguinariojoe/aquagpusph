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
hfac = 3.0
cs = 45.0
courant = 0.25
gamma = 7.0
refd = 998.0
alpha = 0.1
delta = 0.0
visc_dyn = 0.000894
# Tank dimensions
H = 0.508
L = 0.9
D = 0.062
# Fluid
h = 0.093
# Stimated required number of fluid particles
n = 100000

# Dimensions and number of particles readjustment
# ===============================================

nx = int(round((n * D * D / h / L)**(1.0 / 3.0)))
ny = int(round((n * L * L / h / D)**(1.0 / 3.0)))
nz = n // (nx * ny)
n = nx * ny * nz

drx = D/nx
dry = L/ny
drz = h/nz
dr = max(drx, dry, drz)
hFluid = nz*drz
visc_dyn = max(alpha / 10.0 * refd * hfac * dr * cs, visc_dyn)

# Solid boundary elements
Nx = nx
Ny = ny
Nz = int(round(H / drz)) + 1

# Correct tank dimensions
L = Ny * dr
D = Nx * dr
H = Nz * dr

DeLeffeDistFactor = 1
Nx = DeLeffeDistFactor * Nx
Ny = DeLeffeDistFactor * Ny
Nz = DeLeffeDistFactor * Nz
N = 2 * Nx * Ny + 2 * Nx * Nz + 2 * Ny * Nz

# Particles generation
# ====================

prb = cs * cs * refd / gamma

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
#
#    Fluid definition autoexported by AQUAgpusph.
#    Particle needed data are:
#        Position
#        Normal      (Fluid particles can have null normal)
#        Velocity
#        Mass
#    Particle optional data are (sorted):
#        Imove       (1)
#        Density     (refd)
#        Sound speed (cs)
#        KernelH     (h)
#        Ifluid      (ifluid)
#        DensityRate (0.0)
#        Force       ([0,0,0])
#
#    Normal direction is not relevant, but ensure that for
#    solid vertexes are normalised.
#
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
    if Percentage != (i * 100) / n:
        Percentage = (i * 100) / n
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    j = i
    idd = j / ny
    idx = idd / nz
    idy = j % ny
    idz = idd % nz
    imove = 1
    pos = (idx * dr - 0.5 * (D - dr),
           idy * dr - 0.5 * (L - dr),
           idz * dr + 0.5 * dr)
    press = refd * g * (hFluid - pos[2])
    dens = pow(press / prb + 1.0, 1.0 / gamma) * refd
    mass = dens * dr**3.0
    string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
        pos[0], pos[1], pos[2],
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
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
    if Percentage != (i * 100) / N:
        Percentage = (i * 100) / N
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    # Bottom
    if(i < Nx * Ny):
        j = i
        idd = j
        idz = 0
        idx = idd / Ny + 0.5
        idy = idd % Ny + 0.5
        normal = [0.0, 0.0, -1.0]
    # Roof
    elif(i < 2 * Nx * Ny):
        j = i - Nx * Ny
        idd = j
        idz = Nz
        idx = idd / Ny + 0.5
        idy = idd % Ny + 0.5
        normal = [0.0, 0.0, 1.0]
    # Left
    elif(i < 2 * Nx * Ny + Nx * Nz):
        j = i - 2 * Nx * Ny
        idd = j
        idy = 0
        idx = idd / Nz + 0.5
        idz = idd % Nz + 0.5
        normal = [0.0, -1.0, 0.0]
    # Right
    elif(i < 2 * (Nx * Ny + Nx * Nz)):
        j = i - 2 * Nx * Ny - Nx * Nz
        idd = j
        idy = Ny
        idx = idd / Nz + 0.5
        idz = idd % Nz + 0.5
        normal = [0.0, 1.0, 0.0]
    # Front
    elif(i < 2 * (Nx * Ny + Nx * Nz) + Ny * Nz):
        j = i - 2 * (Nx * Ny + Nx * Nz)
        idd = j
        idx = 0
        idy = idd / Nz + 0.5
        idz = idd % Nz + 0.5
        normal = [-1.0, 0.0, 0.0]
    # Back
    elif(i < 2 * (Nx * Ny + Nx * Nz + Ny * Nz)):
        j = i - 2 * (Nx*Ny + Nx * Nz) - Ny * Nz
        idd = j
        idx = nx
        idy = idd / Nz + 0.5
        idz = idd % Nz + 0.5
        normal = [1.0, 0.0, 0.0]
    pos = (idx * dr / DeLeffeDistFactor - 0.5 * D,
           idy * dr / DeLeffeDistFactor - 0.5 * L,
           idz * dr / DeLeffeDistFactor)
    if pos[2] <= hFluid:
        press = refd * g * (hFluid - pos[2])
        dens = pow(press / prb + 1.0, 1.0 / gamma) * refd
    else:
        dens = refd
        press = prb*pow(dens / refd, gamma - 1.0)
    imove = -3
    mass = (dr / DeLeffeDistFactor)**2.0
    string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
        pos[0], pos[1], pos[2],
        normal[0], normal[1], normal[2],
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
print('    100%')
output.close()

# Sensors generation
# ==================

output = open("Sensors.dat", "w")
pos = (0.0, -0.45, 0.093)
normal = (0.0, -1.0, 0.0)
dens = refd
mass = 0.0
imove = 0.0
string = ("{} {} {} 0.0, " * 4 + "{}, {}, {}, {}\n").format(
    pos[0], pos[1], pos[2],
    normal[0], normal[1], normal[2],
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    dens,
    0.0,
    mass,
    imove)
output.write(string)
output.close()

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('Fluids.xml', 'Main.xml', 'Motion.xml', 'Sensors.xml', 'Settings.xml',
       'SPH.xml', 'Time.xml')

radius = math.sqrt((0.5 * L)**2 + H**2)
domain_min = (-D, -1.1 * radius, -0.5 * L * math.sin(0.5 * math.pi), 0.0)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (D, 1.1 * radius, 1.1 * radius, 0.0)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'GAMMA':str(gamma),
        'REFD':str(refd), 'VISC_DYN':str(visc_dyn), 'DELTA':str(delta),
        'G':str(g), 'N_SENSORS':str(1), 'N':str(n + N)}
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