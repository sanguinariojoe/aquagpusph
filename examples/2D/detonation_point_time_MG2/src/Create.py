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

#import os.path as path

import os
import sys
script_folder = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_folder, "../../"))
import aqua_example_utils as utils

import math
import numpy as np


n_H2 = 0.4
n_O2 = 0.21
n_N2 = 0.79

M_H2 = 0.002
M_O2 = 0.032
M_N2 = 0.028

sum_n = n_H2 + n_O2 + n_N2
n_H2 = n_H2 / sum_n
n_O2 = n_O2 / sum_n
n_N2 = n_N2 / sum_n
n_H2O= 0.0
 
MM2 = n_H2 * M_H2 + n_O2 * M_O2 + n_N2 * M_N2

y_H2_2 = n_H2 * M_H2 / MM2
y_O2_2 = n_O2 * M_O2 / MM2
y_N2_2 = n_N2 * M_N2 / MM2
y_H2O_2 = 0.
#z = y_H2

#MM1=0.21 * M_O2 + 0.79 * M_N2
#y_H2_1 = 0.4
#y_O2_1 = 0.21 * M_O2 / MM1
#y_N2_1 = 0.79 * M_N2 / MM1
#y_H2O_1 = 0.

MM1 = MM2
y_H2_1 = y_H2_2
y_O2_1 = y_O2_2
y_N2_1 = y_N2_2
y_H2O_1 = y_H2O_2


print("")
print(f"n_H2 = {n_H2}")
print(f"n_O2 = {n_O2}")
print(f"n_N2 = {n_N2}")
print(f"n_H2O = {n_H2O}")
print("")


courant = 0.1

R = .05
R0 = 0.005

hfac = 4.0

n = 300000

gamma=1.4
cv=716.0

#p1 = 30.0e5
#p1 = 1.0 * 8.31 / MM2 * 6000.0
p1 = 10.0e5
p2 = 1.0e5

print("P1 = %f"%(p1,))
print("P2 = %f"%(p2,))

T1 = 300.0
T2 = 300.0

#rho1 = p1 * MM1 /(8.31 * T1)
rho1 = p1 * MM1 /(8.31 * T1)
#rho2 = 1.00001
rho2 = p2 * MM2 /(8.31 * T2)

print("")
print(f"rho1 = {rho1}")
print(f"rho2 = {rho2}")
print("")

c1 = np.sqrt(gamma * p1 / rho1)
c2 = np.sqrt(gamma * p2 / rho2)

print("")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print("")

ssound=max(c1, c2)

e1=p1 / ((gamma - 1.0) * rho1)
e2=p2 / ((gamma - 1.0) * rho2)
#e_all = p2 / ((gamma - 1.0) * rho_all)
print("")
print(f"e1 = {e1}")
print(f"e2 = {e2}")
print("")
# Distance between particles
# ==========================
Vol = np.pi * R**2
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
                  z=0.0, y_H2=0.0, dy_H2dt=0.0, y_O2=0.0, dy_O2dt=0.0, y_N2=0.0, dy_N2dt=0.0, y_H2O=0.0, dy_H2Odt=0.0, Trigger=0,
                  imove=1):
    m = rho * dr**2
    string = ("{} {}, " * 4 + "{}, {}, {}, {}, {}, " 
              +
              "{} " * 3 + "{}, "   
              +
              "{} " * 3 + "{}, " 
              +
              "{}, {}, {}\n").format(
        p[0], p[1],
        n[0], n[1],
        u[0], u[1],
        dudt[0], dudt[1],  # r, normal, u, dudt,
        rho,
        drhodt,  # rho, drhodt,
        e,
        dedt,  # eint, deintdt,
        z,  # z,
        y_H2, y_O2, y_N2, y_H2O,
        dy_H2dt, dy_O2dt, dy_N2dt, dy_H2Odt, 
        Trigger,
        m,  # m,
 # nu, xi,       
        imove)  # imove
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
        r = np.sqrt(x**2 + y**2)
        if r > R:
            y += dr
            continue

#        rho, ener, y_H2, y_O2, y_N2, y_H2O = (rho1, e1, y_H2_1, y_O2_1, y_N2_1, y_H2O_1) if r < R0 else (rho2, e2, y_H2_2, y_O2_2, y_N2_2, y_H2O_2)

        #rho, ener, y_H2, y_O2, y_N2, y_H2O = (rho_all, e_all, y_H2_2, y_O2_2, y_N2_2, y_H2O_2)
        rho, ener, y_H2, y_O2, y_N2, y_H2O, trigger = (rho1, e1, y_H2_2, y_O2_2, y_N2_2, y_H2O_2, 1) if r < R0 else (rho2, e2, y_H2_2, y_O2_2, y_N2_2, y_H2O_2, 0)
        z=y_H2
        
        #rho, ener, y_H2, y_O2, y_N2, y_H2O = (rho1, e1, y_H2_1, y_O2_1, y_N2_1, y_H2O_1) if r < R0 else (rho2, e2, y_H2_2, y_O2_2, y_N2_2, y_H2O_2)
        writeParticle(output, (x,y), rho=rho, e=ener, z=z, y_H2=y_H2, y_O2=y_O2, y_N2=y_N2, y_H2O=y_H2O, Trigger=trigger)
        N += 1

        y += dr
    x += dr

print(f'{N} particles. Volume = {N * dr**2} vs {Vol}')

# XML definition generation
# =========================

#templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
#XML = ('Fluids.xml', 'Main.xml', 'Settings.xml', 'SPH.xml', 'Time.xml',
#       'BC.xml')


R_domain = R + 4.0 * h
domain_min = (-R_domain, -R_domain)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (R_domain, R_domain)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'H':str(h), 'GAMMA':str(gamma), 'COURANT':str(courant),
        'R':str(R), 'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max,
        'N':str(N), 'DT':str(dt), 'CS':str(ssound), 'DT':str(dt), 'CV':str(cv)}

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
