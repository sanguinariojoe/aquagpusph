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

courant = 0.1
h0 = 0.5
l0 = 0.5
R0 = 0.2

n = 50000

p1=2.0e5
p2=1.0e5

rho1=1.00001
rho2=1.00001

gamma=1.44

c1=np.sqrt(gamma*p1/rho1)
c2=np.sqrt(gamma*p2/rho2)

ssound=max(c1,c2)

e1=p1/((gamma-1.0)*rho1)
e2=p2/((gamma-1.0)*rho2)

#print(e2)
# Distance between particles
# ==========================
Vol = 4 * l0 * h0
dv = Vol / n
#dr = (dv/3.1415)**0.5
dr = dv**0.5/rho1

h = 10.0 * dr
#h=0.1
dt = 1.0E-5

print("")
print("H=%f"%(h,))
print("")

# ================
# ================




# Particles generation
# ====================
def writeParticle(output, p, n=(0.0,0.0), u=(0.0,0.0),
                  dudt=(0.0,0.0), rho=0.0, drhodt=0.0, e=0.0, dedt=0.0, pres=0.0,
                  imove=1):

    m=rho*dr**2
    #m=dr**2
    string = ("{} {}, " * 4 + "{}, {}, {:.4f}, {}, {}, {}, {}\n").format(
        np.float32(p[0]), np.float32(p[1]),
        np.float32(n[0]), np.float32(n[1]),
        np.float32(u[0]), np.float32(u[1]),
        np.float32(dudt[0]), np.float32(dudt[1]),
        np.float32(rho),
        np.float32(drhodt),
        np.float32(e),
        np.float32(dedt),
        np.float32(m),
        np.float32(pres),
        np.int16(imove))
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

x=-l0
y=-h0

while x < l0:
    while y < h0:
        if np.sqrt(x**2 + y**2)<0.45:
            rho=0.0
            ener=0.0
            pres=0.0
            if np.sqrt(x**2 + y**2)<R0:
                rho=rho1
                ener=e1
                pres=p1
            else:
                rho=rho2
                ener=e2
                pres=p2

            writeParticle(output=output, p=(x,y), rho=rho, e=ener, pres=pres, imove=1)
            N += 1
        else:
            writeParticle(output=output, p=(x,y), rho=rho2, e=e2, pres=p2, imove=-1)
            N += 1
        y += dr
    x += dr
    y = -h0


print('{} particles'.format(N))

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('Fluids.xml', 'Main.xml', 'Settings.xml', 'SPH.xml', 'Time.xml')

factor=5.0
domain_min = (-factor*l0, -factor*h0)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (factor*l0, factor*h0)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max,
        'N':str(N), 'DT':str(dt), 'H':str(h), 'CS':str(ssound)}
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

































#
