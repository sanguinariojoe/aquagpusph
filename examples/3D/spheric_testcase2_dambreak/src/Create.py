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
cs = 40.0
courant = 0.25
refd = 998.0
alpha = 0.0
delta = 1.0
visc_dyn = 0.000894
# Initial fluid dimensions
h = 0.55
l = 1.228
d = 1.0
# Box dimensions
H_box = 0.161
L_box = 0.161
D_box = 0.403
x_box = -1.248
# Tank dimensions
H = 1.0
L = l + abs(x_box) + 0.744
D = d
# Stimated required number of fluid particles
n = 100000

Vol = l * d * h
dv = Vol / n
dr = dv**(1.0 / 3.0)

nx = int(round(l / dr))
ny = int(round(d / dr))
nz = int(round(h / dr))

hFluid = nz * dr
visc_dyn = max(alpha / 10.0 * refd * hfac * dr * cs, visc_dyn)

Nx_box = int(round(L_box / dr))
Ny_box = int(round(D_box / dr))
Nz_box = int(round(H_box / dr))
x_box = (int(round(x_box / dr - 0.5 * Nx_box)) + 0.5 * Nx_box) * dr

Nx = int(round(L / dr))
Ny = int(round(D / dr))
Nz = int(round(H / dr))

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

x0 = 0.0
y0 = - 0.5 * Ny * dr
z0 = 0.0

print('Fluid particles...')
n = 0
Percentage = -1
for i in range(nx):
    if Percentage != (n * 100) // (nx * ny * nz):
        Percentage = (n * 100) // (nx * ny * nz)
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    x = x0 + (i + 0.5) * dr
    for j in range(ny):
        y = y0 + (j + 0.5) * dr
        for k in range(nz):
            z = z0 + (k + 0.5) * dr
            imove = 1
            press = refd * g * (hFluid - z) * math.cos(
                0.5 * math.pi * (l - x) / l)
            dens = refd + press / cs**2 
            mass = dens * dr**3.0
            string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
                x, y, z,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                dens,
                0.0,
                mass,
                imove)
            output.write(string)
            n += 1

x0 = x_box - 0.5 * Nx_box * dr
y0 = - 0.5 * Ny_box * dr
z0 = 0.0
dens = refd
press = 0.0
mass = dr**2.0
print('Box boundary elements...')
N_box = 0

print('    Top face...')
normal = (0.0, 0.0, 1.0)
tangent = (-1.0, 0.0, 0.0)
z = z0 + Nz_box * dr
for i in range(Nx_box):
    x = x0 + (i + 0.5) * dr
    for j in range(Ny_box):
        y = y0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N_box += 1

print('    Front face...')
normal = (1.0, 0.0, 0.0)
tangent = (0.0, -1.0, 0.0)
x = x0
for i in range(Ny_box):
    y = y0 + (i + 0.5) * dr
    for j in range(Nz_box):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N_box += 1

print('    Back face...')
normal = (-1.0, 0.0, 0.0)
tangent = (0.0, -1.0, 0.0)
x = x0 + Nx_box * dr
for i in range(Ny_box):
    y = y0 + (i + 0.5) * dr
    for j in range(Nz_box):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N_box += 1

print('    Left face...')
normal = (0.0, 1.0, 0.0)
tangent = (1.0, 0.0, 0.0)
y = y0
for i in range(Nx_box):
    x = x0 + (i + 0.5) * dr
    for j in range(Nz_box):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N_box += 1

print('    Right face...')
normal = (0.0, -1.0, 0.0)
tangent = (-1.0, 0.0, 0.0)
y = y0 + Ny_box * dr
for i in range(Nx_box):
    x = x0 + (i + 0.5) * dr
    for j in range(Nz_box):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N_box += 1

print('OK')

x0 = -(Nx - nx) * dr
y0 = - 0.5 * Ny * dr
z0 = 0.0
dens = refd
press = 0.0
mass = dr**2.0
print('Tank boundary elements...')
N = 0

print('    Bottom face...')
xbox_min = x_box - 0.5 * Nx_box * dr
xbox_max = x_box + 0.5 * Nx_box * dr
ybox_min = -0.5 * Ny_box * dr
ybox_max = 0.5 * Ny_box * dr

normal = (0.0, 0.0, -1.0)
tangent = (-1.0, 0.0, 0.0)
z = z0
for i in range(Nx):
    x = x0 + (i + 0.5) * dr
    for j in range(Ny):
        y = y0 + (j + 0.5) * dr
        if x > xbox_min and x < xbox_max:
            if y > ybox_min and y < ybox_max:
                continue
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Top face...')
normal = (0.0, 0.0, 1.0)
tangent = (1.0, 0.0, 0.0)
z = z0 + Nz * dr
for i in range(Nx):
    x = x0 + (i + 0.5) * dr
    for j in range(Ny):
        y = y0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Front face...')
normal = (-1.0, 0.0, 0.0)
tangent = (0.0, -1.0, 0.0)
x = x0
for i in range(Ny):
    y = y0 + (i + 0.5) * dr
    for j in range(Nz):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Back face...')
normal = (1.0, 0.0, 0.0)
tangent = (0.0, 1.0, 0.0)
x = x0 + Nx * dr
for i in range(Ny):
    y = y0 + (i + 0.5) * dr
    for j in range(Nz):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Left face...')
normal = (0.0, -1.0, 0.0)
tangent = (1.0, 0.0, 0.0)
y = y0
for i in range(Nx):
    x = x0 + (i + 0.5) * dr
    for j in range(Nz):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N += 1

print('    Right face...')
normal = (0.0, 1.0, 0.0)
tangent = (-1.0, 0.0, 0.0)
y = y0 + Ny * dr
for i in range(Nx):
    x = x0 + (i + 0.5) * dr
    for j in range(Nz):
        z = z0 + (j + 0.5) * dr
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            normal[0], normal[1], normal[2],
            tangent[0], tangent[1], tangent[2],
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        N_box += 1

output.close()
print('OK')

positions = []
normals = []
tangents = []

n_sensors = 8
xx = xbox_max
yy = 0.0
zz = Nz_box * dr
nn = (-1.0, 0.0, 0.0)
tt = (0.0, 1.0, 0.0)
for i in range(4):
    normals.append(nn)
    tangents.append(tt)
    positions.append((xx, yy, 0.021 + i * 0.04))

nn = (0.0, 0.0, 1.0)
tt = (1.0, 0.0, 0.0)
for i in range(4):
    normals.append(nn)
    tangents.append(tt)
    positions.append((xx - 0.021 - i * 0.04, yy, zz))

output = open("Sensors.dat", "w")
for i in range(len(positions)):
    pos = positions[i]
    normal = normals[i]
    tangent = tangents[i]
    dens = refd
    mass = 0.0
    imove = 0
    string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
        pos[0], pos[1], pos[2],
        normal[0], normal[1], normal[2],
        tangent[0], tangent[1], tangent[2],
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        dens,
        0.0,
        mass,
        imove)
    output.write(string)
output.close()

s = 10.0 * hfac * dr
domain_min = (-(L - l + s), -(0.5 * D + s), -s, 0.0)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (l + s, 0.5 * D + s, H + s, 0.0)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'N_SENSORS':str(n_sensors), 'N':str(n + N_box + N)}
utils.configure(data, os.path.join(script_folder, "templates"))
