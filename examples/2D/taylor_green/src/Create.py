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

hfac = 4.0
cs = 50.0
courant = 0.1
refd = 1.0
alpha = 0.0
Re = 2000
U = 1.0
L = H = 2.0 * math.pi
periods_to_run = 30.0
p0 = 3.0 * refd * U**2
# Number of fluid particles in y direction
nx = ny = 400

dr = L / nx
n = nx * ny

visc_dyn = refd * U * L / Re
visc_dyn = max(alpha / 8.0 * refd * hfac * dr * cs, visc_dyn)

T = L / U
end_time = periods_to_run * T

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
i = 0
imove = 1
Ekin = 0.0
x = -0.5 * L + 0.5 * dr
while x < 0.5 * L:
    y = -0.5 * H + 0.5 * dr
    while y < 0.5 * H:
        if Percentage != (i * 100) / n:
            Percentage = (i * 100) / n
            if not Percentage % 10:
                string = '    {}%'.format(Percentage)
                print(string)
        
        
        i += 1
        u = U * math.sin(x) * math.cos(y) 
        v = -U * math.cos(x) * math.sin(y)
        press = (math.cos(2.0 * x) + math.cos(2.0 * y)) / 4.0  
        dens = refd + press / cs**2
        mass = refd * dr**2.0
        Ekin += 0.5 * mass * (u*u + v*v)
        string = ("{} {}, " * 4 + "{}, {}, {}, {}\n").format(
            x, y,
            0.0, 0.0,
            u, v,
            0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        y += dr
    x += dr
print('    100%')

domain_min = (-L, -H)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (L, H)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'P0':str(p0), 'RE':str(Re), 'L':str(L),
        'U':str(U), 'N':str(n), 'END_TIME':str(end_time), 'E_KIN':str(Ekin)}
utils.configure(data, os.path.join(script_folder, "templates"))
