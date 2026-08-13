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
import numpy as np


# Tank dimensions
H = 2.0
L = 12.0  # Without the beach
slope_beach = 1 / 15
L_beach = H / slope_beach
# Fluid
h = 1.67
# Wave height probe position
x_probe = 4.0
# Scale factor (in this case it is fixed by the model/real depth ratio)
scale = 1.0 / 18.0
# JONSWAP spectrum
Hs = 3.0
Tp = 10.0

# SPH params
ny = 100  # Number of fluid particles (in z dir)
g = 9.81
hfac = 2.0
cs = 50.0
courant = 0.5
refd = 998.0
alpha = 0.0
delta = 0.1
visc_dyn = 0.000894

# Let's compute the particles interspace and readjust some dimensions
dr = h / ny
Nx = int(round(L / dr))
Nx_beach = int(round(L_beach / dr))
Nx += Nx_beach
nx = Nx
Ny = int(round(H / dr))
L = Nx * dr
L_beach = Nx_beach * dr
H = Ny * dr

length_beach = np.linalg.norm((L_beach, H))
N_beach = int(round(length_beach / dr))
dr_beach = length_beach / N_beach

# Artificial viscosity (if queried)
visc_dyn = max(alpha / 8.0 * refd * hfac * dr * cs, visc_dyn)

print("Writing fluid...")
output = open("Fluid.dat", "w")

# Square flume
x = np.linspace(0.5 * dr, L - 0.5 * dr, num=nx)
y = np.linspace(-h + 0.5 * dr, - 0.5 * dr, num=ny)
xv, yv = np.meshgrid(x, y, indexing='ij')
xv = xv.flatten()
yv = yv.flatten()
# Remove the points on the beach
xb = xv - (L - L_beach)
yb = -h + H * xb / L_beach
mask = np.logical_not(np.logical_and(xb > 0.0, yv < yb + 0.25 * dr))
xv = xv[mask]
yv = yv[mask]
n = xv.shape[0]

Percentage = -1
for i in range(0, n):
    if Percentage != (i * 100) // n:
        Percentage = (i * 100) // n
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    imove = 1
    x, y = xv[i], yv[i]
    press = - refd * g * y
    dens = refd + press / cs**2 
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
print('    100%')
output.close()

output = open("Tank.dat", "w")
string = """
    Writing the solid tank...
"""
print(string)
x = np.linspace(0.5 * dr, L - L_beach - 0.5 * dr, num=Nx - Nx_beach)
y = -np.ones(len(x)) * h
normalx = np.zeros(len(x))
normaly = -np.ones(len(x))
s = np.ones(len(x)) * dr
x_beach = np.linspace(L - L_beach + 0.5 * dr, L - 0.5 * dr, num=N_beach)
y_beach = -h + H * (x_beach - (L - L_beach)) / L_beach
normalx_beach = np.ones(len(x_beach)) * H / length_beach
normaly_beach = -np.ones(len(x_beach)) * L_beach / length_beach
s_beach = np.ones(len(x_beach)) * dr_beach

x = np.concatenate((x, x_beach))
y = np.concatenate((y, y_beach))
normalx = np.concatenate((normalx, normalx_beach))
normaly = np.concatenate((normaly, normaly_beach))
s = np.concatenate((s, s_beach))
n_tank = len(x)

Percentage = -1
for i in range(0, n_tank):
    if Percentage != (i * 100) // n_tank:
        Percentage = (i * 100) // n_tank
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    if y[i] <= 0:
        press = -refd * g * y[i]
        dens = refd + press / cs**2 
    else:
        dens = refd
        press = 0.0
    imove = -3
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        x[i], y[i],
        normalx[i], normaly[i],
        0.0, 0.0,
        0.0, 0.0,
        dens,
        0.0,
        s[i],
        imove)
    output.write(string)
print('    100%')
output.close()

output = open("Flap.dat", "w")
string = """
    Writing the flap...
"""
print(string)
y = np.linspace(-h + 0.5 * dr, H - h - 0.5 * dr, num=Ny)
x = np.zeros(len(y))
n_flap = len(x)

Percentage = -1
for i in range(0, n_flap):
    if Percentage != (i * 100) // n_flap:
        Percentage = (i * 100) // n_flap
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    if y[i] <= 0:
        press = -refd * g * y[i]
        dens = refd + press / cs**2 
    else:
        dens = refd
        press = 0.0
    imove = -3
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        x[i], y[i],
        -1.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        dens,
        0.0,
        dr,
        imove)
    output.write(string)
print('    100%')
output.close()

support = 2 * hfac * dr
domain_min = (-2.0 - 2.0 * support, -h - 2.0 * support)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (L + 2.0 * support, H - h + 2.0 * support)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'L':str(L), 'H':str(H), 'DEPTH':str(h), 'X_PROBE':str(x_probe),
        'SCALE':str(scale), 'Hs':str(Hs), 'Tp':str(Tp),
        'n_fluid':str(n), 'n_tank':str(n_tank), 'n_flap':str(n_flap)}
utils.configure(data, os.path.join(script_folder, "templates"))
