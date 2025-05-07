#! /usr/bin/env python3
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


courant = 0.25
R = 0.5
R0 = 0.2
hfac = 2.0

n = 50000

gamma = 1.4

p1 = 2.0e5
p2 = 1.0e5

rho1 = 1.00001
rho2 = 1.00001

c1 = math.sqrt(gamma * p1 / rho1)
c2 = math.sqrt(gamma * p2 / rho2)

cs = max(c1, c2)

e1 = p1 / ((gamma - 1.0) * rho1)
e2 = p2 / ((gamma - 1.0) * rho2)

t_max = R0 / cs

Vol = math.pi * R**2
dv = Vol / n
dr = dv**0.5
h = hfac * dr

print("")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print(f"dr = {dr}")
print(f"h = {h}")
print("")

def writeParticle(output, p, n=(0.0, 0.0), u=(0.0, 0.0),
                  dudt=(0.0, 0.0), rho=0.0, drhodt=0.0, e=0.0, dedt=0.0,
                  imove=1):
    m = rho * dr**2
    string = ("{} {}, " * 4 + "{}, {}, {}, {}, {}, {}\n").format(
        p[0], p[1],
        n[0], n[1],
        u[0], u[1],
        dudt[0], dudt[1],
        rho,
        drhodt,
        e,
        dedt,
        m,
        imove)
    output.write(string)


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

N = 0
string = """
    Writing reservoir fluid particles...
"""
print(string)

x = -R
while x < R:
    y = -R
    while y < R:
        r = math.sqrt(x**2 + y**2)
        if r > R:
            y += dr
            continue

        rho, ener = (rho1, e1) if r < R0 else (rho2, e2)

        writeParticle(output, (x, y), rho=rho, e=ener)
        N += 1

        y += dr
    x += dr

print(f'{N} particles.')

factor = 2.0
R_domain = R + 4.0 * h
domain_min = (-R_domain, -R_domain)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (R_domain, R_domain)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR': str(dr), 'HFAC': str(hfac), 'H': str(h), 'COURANT': str(courant),
        'R': str(R), 'R0': str(R0), 'T': str(t_max),
        'DOMAIN_MIN': domain_min, 'DOMAIN_MAX': domain_max,
        'N': str(N), 'CS': str(cs), 'GAMMA': str(gamma),
        'P1': str(p1), 'P2': str(p2), 'RHO1': str(rho1), 'RHO2': str(rho2),
        'E1': str(e1), 'E2': str(e2), }
utils.configure(data, os.path.join(script_folder, "templates"))
