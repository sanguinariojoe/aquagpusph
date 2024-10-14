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
hfac = 2.0
cs = 40.0
courant = 0.1
refd = 998.0
alpha = 0.0
delta = 10.0
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

# Parse the available devices
cl_platforms = []
cl_devices = []
cl_types = []
with open("clinfo_platforms", "r") as f1, open("clinfo_devices", "r") as f2:
    for i,line in enumerate(f1.readlines()):
        line = line.strip()
        if line == "":
            continue
        ndevices = int(line.split(" ")[-1])
        for j in range(ndevices):
            cl_platforms.append(str(i))
            cl_devices.append(str(j))
    for line in f2.readlines():
        line = line.strip()
        if line == "":
            continue
        cl_types.append(line.split(" ")[-1])

# Let's try to use just GPUs
if "GPU" in cl_types:
    for i in list(range(len(cl_types)))[::-1]:
        if cl_types[i] != "GPU":
            del cl_platforms[i]
            del cl_devices[i]
            del cl_types[i]

nprocs = len(cl_platforms)
if nprocs == 1:
    # Artificially duplicate the very only available device
    print("WARNING: The same platform and device will be considered twice")
    cl_platforms.append(cl_platforms[0])
    cl_devices.append(cl_devices[0])
    cl_types.append(cl_types[0])
    nprocs = 2

# Create the MPI hostfile
with open("hostfile", "w") as f:
    f.write("localhost slots={0} max-slots={0}\n".format(nprocs))

# Dimensions and number of particles readjustment
# ===============================================

Vol = l * d * h
dv = Vol / n
dr = dv**(1.0 / 3.0)

nx = int(round(l / dr))
ny = int(round(d / dr))
nz = int(round(h / dr))

# Correction to avoid particles laying on processes interface
while (ny - 1) % nprocs == 0:
    ny += 1
    dr = d / ny
nx = int(round(l / dr))
ny = int(round(d / dr))
nz = int(round(h / dr))

hFluid = nz * dr
visc_dyn = max(alpha / 10.0 * refd * hfac * dr * cs, visc_dyn)

Nx_box = int(round(L_box / dr))
Ny_box = int(round(D_box / dr))
Nz_box = int(round(H_box / dr))
x_box = (int(round(x_box / dr - 0.5 * Nx_box)) + 0.5 * Nx_box) * dr

xbox_min = x_box - 0.5 * Nx_box * dr
xbox_max = x_box + 0.5 * Nx_box * dr
ybox_min = -0.5 * Ny_box * dr
ybox_max = 0.5 * Ny_box * dr

Nx = int(round(L / dr))
Ny = int(round(D / dr))
Nz = int(round(H / dr))

domain_min = (-(L - l + dr), -(0.5 * D + dr), -dr, 0.0)
domain_max = (l + dr, 0.5 * D + dr, H + dr, 0.0)

STRING = """#############################################################
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
print(STRING)
def fluid_file(proc):
    output = open("Fluid_{}.dat".format(proc), "w")
    output.write(STRING)

    # Particles generation
    # ====================

    x0 = 0.0
    y0 = - 0.5 * Ny * dr
    z0 = 0.0
    dy = (domain_max[1] - domain_min[1]) / nprocs
    y_min = domain_min[1] + proc * dy
    y_max = domain_min[1] + (proc + 1) * dy

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
            assert (y != y_min) and (y != y_max), \
                "Particle laying in processes interface detected!"
            if not y_min < y < y_max:
                continue
            for k in range(nz):
                z = z0 + (k + 0.5) * dr
                imove = 1
                press = refd * g * (hFluid - z)
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

    # Box generation
    # ==============

    x0 = x_box - 0.5 * Nx_box * dr
    y0 = - 0.5 * Ny_box * dr
    z0 = 0.0
    dens = refd
    press = 0.0
    mass = dr**2.0
    print('Box boundary elements...')
    N_box = 0

    print('    Top face...')
    normal = (0.0, 0.0, -1.0)
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

    # Tank generation
    # ===============

    x0 = -(Nx - nx) * dr
    y0 = - 0.5 * Ny * dr
    z0 = 0.0
    dens = refd
    press = 0.0
    mass = dr**2.0
    print('Tank boundary elements...')
    N = 0

    print('    Bottom face...')
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
            N += 1

    # Buffer particles
    # ================

    print('Buffer particles...')
    x, y, z, _ = domain_max
    press = 0
    dens = refd 
    mass = dens * dr**3.0
    imove = -255
    for i in range(n):
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
    n_buffer = n

    output.close()
    print('OK')

    return n, N_box, N, n_buffer
print('OK')

n = 0
for proc in range(nprocs):
    n += sum(fluid_file(proc))

# Sensors generation
# ==================

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

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')

# Settings is a bit special
with open(path.join(templates_path, 'Settings.xml'), 'r') as fin, \
     open('Settings.xml', 'w') as fout:
    txt = fin.read()
    dev_start = txt.find("        <Device")
    dev_end = dev_start + txt[dev_start:].find("/>") + 3
    dev_str = txt[dev_start:dev_end + 1]
    fout.write(txt[:dev_start])
    for proc in range(nprocs):
        fout.write(dev_str.replace(
            "{{CLPLATFORM}}", cl_platforms[proc]).replace(
            "{{CLDEVICE}}", cl_devices[proc]).replace(
            "{{CLTYPE}}", "ALL"))
    fout.write(txt[dev_end:])

XML = ('Fluids.xml', 'Main.xml', 'MPI.xml', 'Sensors.xml', 'SPH.xml',
       'Time.xml', 'h_sensor.cl')

domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'N_SENSORS':str(n_sensors), 'N':str(n)}
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
