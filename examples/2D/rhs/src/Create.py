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
import numpy as np

#zone 2 h2+air
nu=1.0e-5
xi=1.0e-2

n_H2 = 0.0
n_O2 = 0.21
n_N2 = 0.79

M_H2 = 0.002
M_O2 = 0.032
M_N2 = 0.028

sum_n = n_H2 + n_O2 + n_N2
n_H2 = n_H2 / sum_n
n_O2 = n_O2 / sum_n
n_N2 = n_N2 / sum_n

MM2 = n_H2 * M_H2 + n_O2 * M_O2 + n_N2 * M_N2

y_H2_2 = n_H2 * M_H2 / MM2
y_O2_2 = n_O2 * M_O2 / MM2
y_N2_2 = n_N2 * M_N2 / MM2
y_H2O_2 = 0.
#z = y_H2

#zone 1 pure air

MM1=0.21 * M_O2 + 0.79 * M_N2
y_H2_1 = 0.0
y_O2_1 = 0.21 * M_O2 / MM1
y_N2_1 = 0.79 * M_N2 / MM1
y_H2O_1 = 0.

courant = 0.1
support = 2.0 

b0plus=0.1
b0minus=0.1

hfac = 2.0

nx = 100

p1 = 1.0e5
p2 = 1.0e5

T1 = 2000.0

rho1 = p1 * MM1 /(8.31 * T1)
rho2 = 1.00001

gamma=1.4

c1 = np.sqrt(gamma * p1 / rho1)
c2 = np.sqrt(gamma * p2 / rho2)

print("")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print("")

ssound=max(c1, c2)
cs = ssound
e1=p1 / ((gamma - 1.0) * rho1)
e2=p2 / ((gamma - 1.0) * rho2)

# Distance between particles
# ==========================
courant = 0.5
L = b0minus+b0plus
dr = L / nx
h = hfac * dr
R = 2 * support * h + dr
#dt = 1.0E-5
dt = 0.1*min(dr / c1 , dr / c2)
t_max = (0.5 * L - support * h) / cs
t_max = 1.0e-3

print("")
print(f"dr = {dr}")
print(f"h = {h}")
print(f"R={R}")
print("")


# Particles generation
# ====================
def writeParticle(output, p, n=(0.0, 0.0), u=(0.0, 0.0),
                  dudt=(0.0, 0.0), rho=0.0, drhodt=0.0, e=0.0, dedt=0.0,
                  z=0, y_H2=0, dy_H2dt=0, y_O2=0, dy_O2dt=0, y_N2=0, dy_N2dt=0, y_H2O=0, dy_H2Odt=0,
                  imove=1):
    m = rho * dr**2
    string = ("{} {}, " * 4 + "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n").format(
        p[0], p[1],
        n[0], n[1],
        u[0], u[1],
        dudt[0], dudt[1],#r, normal, u, dudt, 
        rho,
        drhodt,#rho, drhodt,
        e,
        dedt,#eint, deintdt, 
        z,#z, 
        y_H2, dy_H2dt,#y_H2, dy_H2dt, 
        y_O2, dy_O2dt,#y_O2, dy_O2dt, 
        y_N2, dy_N2dt,#y_N2, dy_N2dt, 
        y_H2O, dy_H2Odt,#y_H2O, dy_H2Odt, 
        m,#m, 
        nu,
        xi,#nu, xi, 
        imove)#imove
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

x = -b0minus + 0.5 * dr 
while x < b0plus:
    y = -0.5 * (R - dr)
    while y < 0.5 * R:

        rho, ener, y_H2, y_O2, y_N2, y_H2O = (rho1, e1, y_H2_1, y_O2_1, y_N2_1, y_H2O_1) if x < 0 else (rho2, e2, y_H2_2, y_O2_2, y_N2_2, y_H2O_2)
        z=y_H2
        
        writeParticle(output, (x,y), rho=rho, e=ener, z=z, y_H2=y_H2, y_O2=y_O2, y_N2=y_N2, y_H2O=y_H2O)
        N += 1
        y += dr
    x += dr

#print(f'{N} particles. Volume = {N * dr**2} vs {Vol}')

print(f'{N} particles. Volume = {N * dr**2}')
# Domain definition & Symmetry particles
# ==================

domain_min = (-b0minus - support * h, -1.5 * R - support * h)
domain_max = (b0plus + support * h, 1.5 * R + support * h)

x, y = domain_max[0] + support * h, domain_max[1] + support * h

for _ in range(2*N):
    rho, ener, y_H2, y_O2, y_N2, y_H2O = (rho1, e1, y_H2_1, y_O2_1, y_N2_1, y_H2O_1)
    z=y_H2
    writeParticle(output, (x,y), rho=rho, e=ener, z=z, y_H2=y_H2, y_O2=y_O2, y_N2=y_N2, y_H2O=y_H2O, imove=-255)
N *= 3

domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = str(domain_max).replace('(', '').replace(')', '')

# XML definition generation
# =========================

#templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
#XML = ('Fluids.xml', 'Main.xml', 'Settings.xml', 'SPH.xml', 'Time.xml',
#       'BC.xml')

#data = {'DR':str(dr), 'HFAC':str(hfac), 'H':str(h),  'COURANT':str(courant),
#        'B0minus':str(b0minus),'B0plus':str(b0plus), 'R': str(R), 
#        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max,
#        'N':str(N),  'CS':str(ssound), 'DT':str(dt)}

data = {'DR': str(dr), 'HFAC': str(hfac), 'H': str(h), 'COURANT': str(courant),
        'L': str(L), 'R': str(R), 'T': str(t_max),
        'DOMAIN_MIN': domain_min, 'DOMAIN_MAX': domain_max,
        'N': str(N), 'CS': str(cs), 'GAMMA': str(gamma), 
        'P1': str(p1), 'P2': str(p2), 'RHO1': str(rho1), 'RHO2': str(rho2),
        'E1': str(e1), 'E2': str(e2), }
utils.configure(data, os.path.join(script_folder, "templates"))

#for fname in XML:
#    # Read the template
#    f = open(path.join(templates_path, fname), 'r')
#    txt = f.read()
#    f.close()
#    # Replace the data
#    for k in data.keys():
#        txt = txt.replace('{{' + k + '}}', data[k])
#    # Write the file
#    f = open(fname, 'w')
#    f.write(txt)
#    f.close()
