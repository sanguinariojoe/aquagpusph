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
courant = 0.2
cs_max = 500.0
young = 0.05e6
poisson = 0.49
refd = 2.0e-3
time_to_run = 10.0
# Fluid dimensions
L = 1.0
H = 0.1 * L
# Number of solid particles in y direction
ny = 10
nx = 10 * ny

# Sound speed and shear modulus
nu = poisson
E = young
mu = E / (2.0 * (1.0 + nu))
cs = math.sqrt(mu / refd * (2.0 + 4.0 * nu) / (3.0 - 6.0 * nu))
if cs > cs_max:
    cs = cs_max
    K = refd * cs**2
    mu = 3.0 * K * (1.0 - 2.0 * nu) / (1.0 + 2.0 * nu) 
    E = mu * 2.0 * (1.0 + nu)

# Dimensions and number of particles readjustment
# ===============================================
dr = H / ny
n = nx * ny

end_time = time_to_run

# Particles generation
# ====================
input_template = "{} {}, " * 4 + "{}, " * 2 + "{} {} {} {}," * 2 + "{}, {}\n"

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

print("")
print("E = {}".format(E))
print("nu = {}".format(nu))
print("mu = {}".format(mu))
print("cs = {}".format(cs))
print("")

print("Writing volume particles...")
n_vol = 0
output = open("Solid.dat", "w")
output.write(header)
Percentage = -1
for i in range(nx):
    for j in range(ny):
        id = j + i * ny
        if Percentage != (id * 100) / n:
            Percentage = (id * 100) / n
            if not Percentage % 10:
                string = '    {}%'.format(Percentage)
                print(string)
        imove = 2
        pos = (i * dr + 0.5 * dr,
               j * dr + 0.5 * dr - 0.5 * H)
        press = refd * g * (0.5 * H - pos[1]) 
        dens = refd + press / cs**2
        mass = refd * dr**2.0
        string = input_template.format(
            pos[0], pos[1],      # r
            0.0, 0.0,            # normal
            0.0, 0.0,            # u
            0.0, 0.0,            # dudt
            dens,                # rho
            0.0,                 # drhodt
            0.0, 0.0, 0.0, 0.0,  # S
            0.0, 0.0, 0.0, 0.0,  # dSdt
            mass,                # m
            imove)               # imove
        output.write(string)
        n_vol += 1
print('    100%')
output.close()

print("Writing built in boundary...")
n_left = 0
output = open("BuiltInBC.dat", "w")
output.write(header)
Percentage = -1
for i in range(ny):
    if Percentage != (i * 100) / ny:
        Percentage = (i * 100) / ny
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    imove = -3
    pos = (0.0,
           i * dr + 0.5 * dr - 0.5 * H)
    normal = (-1.0, 0.0)
    press = refd * g * (0.5 * H - pos[1]) 
    dens = refd + press / cs**2
    mass = dr
    string = input_template.format(
        pos[0], pos[1],        # r
        normal[0], normal[1],  # normal
        0.0, 0.0,              # u
        0.0, 0.0,              # dudt
        dens,                  # rho
        0.0,                   # drhodt
        0.0, 0.0, 0.0, 0.0,    # S
        0.0, 0.0, 0.0, 0.0,    # dSdt
        mass,                  # m
        imove)                 # imove
    output.write(string)
    n_left += 1
print('    100%')
output.close()

print("Writing free boundary...")
n_free = 0
output = open("FreeBC.dat", "w")
output.write(header)
Percentage = -1
for i in range(ny):
    if Percentage != (n_free * 100) / (ny + 2 * nx):
        Percentage = (n_free * 100) / (ny + 2 * nx)
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    imove = -3
    pos = (L,
           i * dr + 0.5 * dr - 0.5 * H)
    normal = (1.0, 0.0)
    press = refd * g * (0.5 * H - pos[1]) 
    dens = refd + press / cs**2
    mass = dr
    string = input_template.format(
        pos[0], pos[1],        # r
        normal[0], normal[1],  # normal
        0.0, 0.0,              # u
        0.0, 0.0,              # dudt
        dens,                  # rho
        0.0,                   # drhodt
        0.0, 0.0, 0.0, 0.0,    # S
        0.0, 0.0, 0.0, 0.0,    # dSdt
        mass,                  # m
        imove)                 # imove
    output.write(string)
    n_free += 1
for i in range(nx):
    if Percentage != (n_free * 100) / (ny + 2 * nx):
        Percentage = (n_free * 100) / (ny + 2 * nx)
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    imove = -3
    # Bottom particle
    pos = [i * dr + 0.5 * dr,
           -0.5 * H]
    normal = [0.0, -1.0]
    press = refd * g * (0.5 * H - pos[1]) 
    dens = refd + press / cs**2
    mass = dr
    string = input_template.format(
        pos[0], pos[1],        # r
        normal[0], normal[1],  # normal
        0.0, 0.0,              # u
        0.0, 0.0,              # dudt
        dens,                  # rho
        0.0,                   # drhodt
        0.0, 0.0, 0.0, 0.0,    # S
        0.0, 0.0, 0.0, 0.0,    # dSdt
        mass,                  # m
        imove)                 # imove
    output.write(string)
    n_free += 1
    #Top particle
    pos[1] *= -1.0
    normal[1] *= -1.0
    press = refd * g * (0.5 * H - pos[1]) 
    dens = refd + press / cs**2
    string = input_template.format(
        pos[0], pos[1],        # r
        normal[0], normal[1],  # normal
        0.0, 0.0,              # u
        0.0, 0.0,              # dudt
        dens,                  # rho
        0.0,                   # drhodt
        0.0, 0.0, 0.0, 0.0,    # S
        0.0, 0.0, 0.0, 0.0,    # dSdt
        mass,                  # m
        imove)                 # imove
    output.write(string)
    n_free += 1
print('    100%')
output.close()

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('BCs.xml', 'Main.xml', 'Settings.xml', 'Solid.xml', 'SPH.xml',
       'Time.xml')

domain_min = (-0.25 * L, -H)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (1.2 * L, H)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd), 
        'MU':str(mu), 'G':str(g), 'NSOLID':str(n_vol),
        'NBUILTIN':str(n_left), 'NFREE':str(n_free),
        'L':str(L), 'H':str(H), 'END_TIME':str(end_time)}
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
