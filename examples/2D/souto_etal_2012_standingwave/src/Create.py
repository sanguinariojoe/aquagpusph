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

g = 1
hfac = 4.0
cs = 50.0
courant = 0.1
refd = 1.0
alpha = 0.0
delta = 10.0
periods_to_run = 8
# Fluid dimensions
L = 2.0
H = 0.5 * L
# Wave
Lambda = L
k = 2.0 * math.pi / Lambda
omega = math.sqrt(g * k * math.tanh(k * H))
Epsilon = 0.1
Re = 250
# Number of fluid particles in y direction
ny = 100
nx = 2 * ny

sep = 2.0
dr = H / ny
h = hfac * dr
n = nx * ny

visc_dyn = refd * H * math.sqrt(g * H) / Re
visc_dyn = max(alpha / 8.0 * refd * h * cs, visc_dyn)
A = 0.5 * Epsilon * H

domain_min = (-0.2 * L - 6.0 * sep * h, -0.2 * (H + A) - 6.0 * sep * h)
domain_max = (1.2 * L + 6.0 * sep * h, 1.2 * (H + A) + 6.0 * sep * h)

# Solid boundary elements
Nx = nx

DeLeffeDistFactor = 1
Nx = DeLeffeDistFactor * Nx
N = Nx

# Buffer particles for the symmetry planes
n_buffer = n + N

T = 2.0 * math.pi / omega
end_time = periods_to_run * T

t = 0.0
dt = T / 16
t_list = []
t_adim_list = []
nu = H * math.sqrt(g * H) / Re
E_list = []
E_adim_list = []
while t <= end_time:
    t_list.append(t)
    t_adim_list.append(t / T)
    E_list.append(Epsilon**2 * g * H**2 * Lambda / 32 *
                  # (1 + math.cos(2 * omega * t)) *
                  2 *
                  math.exp(-4 * nu * k**2 * t))
    E_adim_list.append(E_list[-1] / E_list[0])
    t += dt
output = open("Ekin.dat", "w")
output.write("# t\tt/T\tEkin(t)\tEkin(t)/Ekin(0)")
for i in range(len(t_list)):
    output.write("{}\t{}\t{}\t{}\n".format(t_list[i],
                                           t_adim_list[i],
                                           E_list[i],
                                           E_adim_list[i]))
output.close()
Ekin0 = E_list[0]
Epot0 = 0.5 * g * H * (L * H * refd)

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
    pos = (idx * dr + 0.5 * dr,
           idy * dr + 0.5 * dr)
    x = pos[0]
    y = pos[1] - H
    k_u = Epsilon * g * H * k / (2.0 * omega * math.cosh(k * H))
    vel = (k_u * math.sin(k * x) * math.cosh(k * (H + y)),
          -k_u * math.cos(k * x) * math.sinh(k * (H + y)))
    press = refd * g * (H - pos[1])
    dens = refd + press / cs**2 
    mass = dens * dr**2.0
    string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
        pos[0], pos[1],
        0.0, 0.0,
        vel[0], vel[1],
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
    j = i
    idx = j + 0.5
    normal = [0.0, -1.0]
    pos = (idx * dr / DeLeffeDistFactor, 0.0)
    press = refd * g * (H - pos[1])
    dens = refd + press / cs**2 
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
    imove = -255
    dens = refd
    mass = 0.0
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
output.close()

domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'N':str(n + N + n_buffer), 'L':str(L), 'END_TIME':str(end_time),
        'E_KIN':str(Ekin0), 'E_POT':str(Epot0)}
utils.configure(data, os.path.join(script_folder, "templates"))
