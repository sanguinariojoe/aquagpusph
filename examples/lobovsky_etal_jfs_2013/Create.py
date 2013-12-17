#! /usr/bin/env python
#########################################################
#                                                       #
#    #    ##   #  #   #                           #     #
#   # #  #  #  #  #  # #                          #     #
#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #
#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #
#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #
#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #
#                            # #             #          #
#                          ##  #             #          #
#                                                       #
#########################################################
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
#########################################################
#
# Create.py
#
# Script from example lobovsky_etal_jfs_2013. This script can create
# the fluid and the solid file for this case, in .dat format.
#
# This case is a dam break problem with a good experimental data, which
# can be found at arxiv:
#     http://arxiv.org/pdf/1308.0115v1.pdf
# In this example the case of 300mm of initial water level is selected,
# and the gate removing process is not simulated.
#
# Input data:
#
#     h  = Fluid height
#     l  = Fluid length
#     H  = Tank height
#     L  = Tank length
#     n  = Stimated number of particles
#     cs = Sound speed
#     refd  = Density reference
#     gamma = Exponent for Batchelor'67 state equation
#
#########################################################

import math

# General settings
cs     = 45.0
refd   = 997.0
gamma  = 7.0
g      = 9.81
# Tank dimensions
H      = 0.6
L      = 1.61
# Fluid
h      = 0.3
l      = 0.6
# Stimated required number of fluid particles
n      = 500000
# Calculate the number of particles at X/Y
ny     = int(round(math.sqrt(n*h/l)))
nx     = int(round(math.sqrt(n*l/h)))
n      = int(nx*ny)
# Calculate distance between particles
dr     = l/nx
hFluid = (ny-0.5)*dr
# Calculate the number of solid particles
Nx     = int(round(L/dr)) + 1
Ny     = int(round(H/dr)) + 1
Drx    = L / Nx
Dry    = H / Ny
# Correct DeLeffe vertices distance
DeLeffeDistFactor = 1
Nx     = DeLeffeDistFactor*Nx
Ny     = DeLeffeDistFactor*Ny
N      = Nx + 2*Ny
# Some fluid properties
prb    = cs*cs*refd/gamma;
# Open output file, and write initial data
print("Opening output file...")
output = open("Fluid.dat","w")
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
#        Normal      (Fluid particles can have null normal)
#        Velocity
#        Mass
#    Particle optional data are (sorted):
#        Imove       (1)
#        Density     (refd)
#        Sound speed (cs)
#        KernelH     (h)
#        Ifluid      (ifluid)
#        DensityRate (0.0)
#        Force       ([0,0,0])
#
#    Normal direction is not relevant, but ensure that for
#    solid vertexes are normalised.
#
#############################################################
"""
output.write(string)
print(string)
# Write data
string = """
    Writing fluid particles...
"""
print(string)
Percentage=-1
for i in range(0,n):
	if Percentage != (i*100) / n:
		Percentage = (i*100) / n
		if not Percentage%10:
			string = '    {0}%'.format(Percentage)
			print(string)
	j      = i
	idy    = j/nx
	idx    = j%nx
	imove  = 1
	pos    = (0.5*dr + idx*dr, 0.5*dr + idy*dr)
	press  = refd*g*(hFluid-pos[1]);
	dens   = pow( press/prb + 1.0, 1.0/gamma )*refd;
	mass   = dr*dr*dens
	string = """{0} {1} 0.0 0.0 0.0 0.0 {2} {3} {4} {5}\n""".format(pos[0], pos[1], mass, imove, dens, cs)
	output.write(string)
print('    100%')
# Write data
string = """
    Writing nodes...
"""
print(string)
Percentage=-1
for i in range(0,N):
	if Percentage != (i*100) / N:
		Percentage = (i*100) / N
		if not Percentage%10:
			string = '    {0}%'.format(Percentage)
			print(string)
	if(i < Nx):                 # Bottom
		j = i
		idy = 0.0
		idx = j+0.5
		normal = [0.0,-1.0]
		mass  = Drx/DeLeffeDistFactor
	elif(i < Nx + Ny):          # Right
		j = i-Nx
		idx = Nx
		idy = j+0.5
		normal = [1.0,0.0]
		mass  = Dry/DeLeffeDistFactor
	elif(i < Nx + 2*Ny):        # Left
		j = i-Nx-Ny
		idx = 0.0
		idy = j+0.5
		normal = [-1.0,0.0]
		mass  = Dry/DeLeffeDistFactor
	pos = (idx*Drx/DeLeffeDistFactor - L + l, idy*Dry/DeLeffeDistFactor)
	if pos[1] <= h:
		press = refd*g*(hFluid-pos[1]);
		dens  = pow( press/prb + 1.0, 1.0/gamma )*refd;
	else:
		dens  = refd;
		press = prb*pow( dens/refd , gamma - 1.0 );
	imove = -1
	string = """{0} {1} {2} {3} 0.0 0.0 {4} {5} {6} {7}\n""".format(pos[0], pos[1], normal[0], normal[1], mass, imove, dens, cs)
	output.write(string)
print('    100%')
string = """
{0} particles has been written into file ({1} of fluid).
PLEASE, CHANGE THIS VALUE IN THE FLUID XML DEFINITION.
""".format(n+N, n)
print(string)
# Print some useful data
print('Particle distance:  {0}'.format(dr))
print('PLEASE, CHANGE THIS VALUE IN THE SPH XML DEFINITION.')
print('Fluid height error: {0}'.format(abs(h - ny*dr)))
# Clean up
output.close()
