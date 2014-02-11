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

# Generate the fluid particles and the boundary elements for the example
# LateralWater_1x_DeLeffe, Correponding to the lateral water impact SPHERIC
# test case number 10. See
# doc/SOUTOIGLESIAS_BOTIA_SPHERIC_TESTCASE_10.pdf
# In this example boundary integrals are applied for the boundary conditions.

import math

# Input data
# ==========

cs = 45.0
refd = 998.0
gamma = 1.0
g = 9.81
# Tank dimensions
L = 0.9
H = 0.508
# Fluid dimensions
h = 0.093
# Stimated required number of fluid particles
n = 50000

# Dimensions and number of particles readjustment
# ===============================================

ny = int(round(math.sqrt(n * h / L)))
nx = int(round(math.sqrt(n * L / h)))
n = int(nx * ny)

dr = L / nx
hFluid = (ny - 0.5) * dr

# Solid boundary elements
Nx = nx
Ny = int(round(H / dr)) + 1

L = Nx * dr
H = Ny * dr

DeLeffeDistFactor = 1
Nx = DeLeffeDistFactor * Nx
Ny = DeLeffeDistFactor * Ny
N = 2 * Nx + 2 * Ny

# Fluid generation
# ================
prb = cs * cs * refd / gamma

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
#
#    Fluid definition autoexported by AQUAgpusph.
#    Particle needed data are:
#        Position
#        Normal      (Fluid particles can have a null normal)
#        Velocity
#        Mass
#    Particle optional data are (sorted):
#        Imove       (1)
#        Density     (refd)
#        KernelH     (h)
#        Ifluid      (ifluid)
#
#############################################################
"""
output.write(string)
print(string)

string = """
    Writing the fluid particles...
"""
print(string)

Percentage = -1
for i in range(n):
    if Percentage != (i * 100) / n:
        Percentage = (i * 100) / n
        if not Percentage % 10:
            print('    {}%'.format(Percentage))
    j = i
    idy = j / nx
    idx = j % nx
    imove = 1
    pos = (idx * dr - 0.5 * (L - dr), idy * dr + 0.5 * dr)
    press = refd * g * (hFluid - pos[1])
    dens = pow(press / prb + 1.0, 1.0 / gamma) * refd
    mass = dr * dr * dens
    string = "{} {} {} {} {} {} {} {} {} {} {} {}\n".format(pos[0],
                                                            pos[1],
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            dens,
                                                            0.0,
                                                            mass,
                                                            imove)
    output.write(string)
print('    100%%')

string = """
    Writing the boundary elements...
"""
print(string)

Percentage = -1
for i in range(N):
    if Percentage != (i * 100) / N:
        Percentage = (i * 100) / N
        if not Percentage % 10:
            string = '    {}%'.format(Percentage)
            print(string)
    # Bottom face
    if(i < Nx):
        j = i
        idy = 0.0
        idx = j + 0.5
        normal = [0.0, -1.0]
    # Right face
    elif(i < Nx + Ny):
        j = i - Nx
        idx = Nx
        idy = j + 0.5
        normal = [1.0, 0.0]
    # Top face
    elif(i < 2 * Nx + Ny):
        j = i - Nx - Ny
        idy = Ny
        idx = Nx - j - 0.5
        normal = [0.0, 1.0]
    # Left face
    elif(i < 2 * Nx + 2 * Ny):
        j = i - 2 * Nx - Ny
        idx = 0.0
        idy = Ny - j - 0.5
        normal = [-1.0, 0.0]
    pos = (idx * dr / DeLeffeDistFactor - 0.5 * L,
           idy * dr / DeLeffeDistFactor)
    if pos[1] <= h:
        press = refd * g * (hFluid - pos[1])
        dens = pow(press / prb + 1.0, 1.0 / gamma) * refd
    else:
        dens = refd
        press = prb * pow(dens / refd, gamma - 1.0)
    imove = -1
    mass = dr / DeLeffeDistFactor
    string = "{} {} {} {} {} {} {} {} {} {} {} {}\n".format(pos[0],
                                                            pos[1],
                                                            normal[0],
                                                            normal[1],
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            dens,
                                                            0.0,
                                                            mass,
                                                            imove)
    output.write(string)
print('    100%')
string = """
{} particles has been written into file ({} of fluid).
PLEASE, CHANGE THIS VALUE IN THE FLUID XML DEFINITION.
""".format(n + N, n)
print(string)
# Print some useful data
print('Particle distance:  {}'.format(dr))
print('PLEASE, CHANGE THIS VALUE IN THE SPH XML DEFINITION.')
print('Fluid height error: {}'.format(abs(h - ny * dr)))
# Clean up
output.close()
