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

import os.path as path
import math
import numpy as np


n_H2 = 0.2
n_O2 = 0.21
n_N2 = 0.79

M_H2 = 0.002
M_O2 = 0.032
M_N2 = 0.028

sum_n = n_H2 + n_O2 + n_N2
n_H2 = n_H2 / sum_n
n_O2 = n_O2 / sum_n
n_N2 = n_N2 / sum_n

MM1 = n_H2 * M_H2 + n_O2 * M_O2 + n_N2 * M_N2

y_H2_2 = n_H2 * M_H2 / MM1
y_O2_2 = n_O2 * M_O2 / MM1
y_N2_2 = n_N2 * M_N2 / MM1
y_H2O_2 = 0.
#z = y_H2

MM2=0.21 * M_O2 + 0.79 * M_N2
y_H2_1 = 0.0
y_O2_1 = 0.21 * M_O2 / MM2
y_N2_1 = 0.79 * M_N2 / MM2
y_H2O_1 = 0.

courant = 0.1

h0=0.1
b0plus=0.7
b0minus=0.1

hfac = 2.0

n = 50000

gamma=1.4
cv=716.0

p1 = 30.0e5
p2 = 1.0e5

T1 = 2000.0

rho1 = p1 * MM2 /(8.31 * T1)
rho2 = 1.00001

c1 = np.sqrt(gamma * p1 / rho1)
c2 = np.sqrt(gamma * p2 / rho2)

print("")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print("")

ssound=max(c1, c2)

e1=p1 / ((gamma - 1.0) * rho1)
e2=p2 / ((gamma - 1.0) * rho2)

# Distance between particles
# ==========================
Vol = 2.0 * h0 * b0minus + 2.0 * h0 * b0plus 
dv = Vol / n
dr = dv**0.5
h = hfac * dr

#dt = 1.0E-5
dt = 0.1*min(dr / c1 , dr / c2)

print("")
print(f"dr = {dr}")
print(f"h = {h}")
print("")


# Particles generation
# ====================
def writeParticle(output, p, n=(0.0, 0.0), u=(0.0, 0.0),
                  dudt=(0.0, 0.0), rho=0.0, drhodt=0.0, e=0.0, dedt=0.0,
                  z=0, y_H2=0, dy_H2dt=0, y_O2=0, dy_O2dt=0, y_N2=0, dy_N2dt=0, y_H2O=0, dy_H2Odt=0,
                  imove=1):
    m = rho * dr**2
    string = ("{} {}, " * 4 + "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n").format(
        p[0], p[1],
        n[0], n[1],
        u[0], u[1],
        dudt[0], dudt[1],
        rho,
        drhodt,
        e,
        dedt,
        z,
        y_H2, dy_H2dt,
        y_O2, dy_O2dt,
        y_N2, dy_N2dt,
        y_H2O, dy_H2Odt,
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

x = -b0minus
while x < b0plus:
    y = -h0
    while y < h0:

        rho, ener, y_H2, y_O2, y_N2, y_H2O = (rho1, e1, y_H2_1, y_O2_1, y_N2_1, y_H2O_1) if x < 0 else (rho2, e2, y_H2_2, y_O2_2, y_N2_2, y_H2O_2)
        z=y_H2
        writeParticle(output, (x,y), rho=rho, e=ener, z=z, y_H2=y_H2, y_O2=y_O2, y_N2=y_N2, y_H2O=y_H2O)
        N += 1

        y += dr
    x += dr

print(f'{N} particles. Volume = {N * dr**2} vs {Vol}')

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('Fluids.xml', 'Main.xml', 'Settings.xml', 'SPH.xml', 'Time.xml',
       'BC.xml')

H_domain = h0 + 4.0 * h
B_domain_minus = b0minus + 4.0 * h
B_domain_plus = b0plus + 4.0 * h

domain_min = (-B_domain_minus, -H_domain)
domain_min = str(domain_min).replace('(', '').replace(')', '')

domain_max = (B_domain_plus, H_domain)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'H':str(h), 'GAMMA':str(gamma), 'COURANT':str(courant),
        'B0minus':str(b0minus),'B0plus':str(b0plus), 'H0':str(h0), 'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max,
        'N':str(N), 'DT':str(dt), 'CS':str(ssound), 'DT':str(dt), 'CV':str(cv)}

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
