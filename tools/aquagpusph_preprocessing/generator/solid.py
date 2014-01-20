#******************************************************************************
#                                                                             *
#               *    **   *  *   *                           *                *
#              * *  *  *  *  *  * *                          *                *
#             ***** *  *  *  * *****  **  ***  *  *  ** ***  ***              *
#             *   * *  *  *  * *   * *  * *  * *  * *   *  * *  *             *
#             *   * *  *  *  * *   * *  * *  * *  *   * *  * *  *             *
#             *   *  ** *  **  *   *  *** ***   *** **  ***  *  *             *
#                                       * *             *                     *
#                                     **  *             *                     *
#                                                                             *
#******************************************************************************
#                                                                             *
#  This file is part of AQUAgpusph, a free CFD program based on SPH.          *
#  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>               *
#                                                                             *
#  AQUAgpusph is free software: you can redistribute it and/or modify         *
#  it under the terms of the GNU General Public License as published by       *
#  the Free Software Foundation, either version 3 of the License, or          *
#  (at your option) any later version.                                        *
#                                                                             *
#  AQUAgpusph is distributed in the hope that it will be useful,              *
#  but WITHOUT ANY WARRANTY; without even the implied warranty of             *
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
#  GNU General Public License for more details.                               *
#                                                                             *
#  You should have received a copy of the GNU General Public License          *
#  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.        *
#                                                                             *
#******************************************************************************

from numpy import *
from vec import *
from sys import stdout

def perform(faces, level, dr, refd, cs, gamma, g):
	""" Create solid particles
	@param faces Faces defined by 4 points and a normal
	@param level Fluid level
	@param dr Distance between particles
	@param refd Fluid density
	@param cs Sound speed
	@param gamma Batchelor's 67 compressibility factor
	@param g Gravity acceleration
	@return Boundary vertexes
	"""
	parts      = []
	prb        = cs*cs*refd/gamma
	percentage = 0
	idf        = 0
	stdout.write("%d%%..." % (percentage))
	stdout.flush()
	for face in faces:
		# Compute the number of points in each direction
		u0  = face[1] - face[0]
		u1  = face[2] - face[3]
		v0  = face[3] - face[0]
		v1  = face[2] - face[1]
		lu0 = length(u0)
		lu1 = length(u1)
		lv0 = length(v0)
		lv1 = length(v1)
		# nu  = (int(lu0/dr) + int(lu1/dr))/2
		# nv  = (int(lv0/dr) + int(lv1/dr))/2
		nu  = int(max(ceil(lu0/dr), ceil(lu1/dr)))
		nv  = int(max(ceil(lv0/dr), ceil(lv1/dr)))
		du  = 1.0 / nu
		dv  = 1.0 / nv
		# Compute the face are from the triangles area
		ld0 = sqrt(lu0*lu0 + lv0*lv0)
		ld1 = sqrt(lu1*lu1 + lv1*lv1)
		sp0 = 0.5*(lu0 + lv0 + ld0)
		sp1 = 0.5*(lu1 + lv1 + ld1)
		a0  = sqrt(sp0*(sp0 - lu0)*(sp0 - lv0)*(sp0 - ld0))
		a1  = sqrt(sp1*(sp1 - lu1)*(sp1 - lv1)*(sp1 - ld1))
		A   = a0 + a1
		# Create the area elements
		a   = 0.0
		for i in range(0,nu):
			for j in range(0,nv):
				# Get point
				ufact = (i+0.5)/nu
				vfact = (j+0.5)/nv
				point = face[0] + ufact*u0 + vfact*v0 + ufact*vfact*(u1 - u0)
				# Compute area associated to this point
				dpdu  = u0 + vfact*(u1-u0)
				dpdv  = v0 + ufact*(u1-u0)
				area  = linalg.norm( cross(dpdu,dpdv) ) * du * dv
				a     = a + area
				# Compute density from hidrostatic
				if point[2] <= level:
					press = refd*g*(level-point[2])
					dens  = pow( press/prb + 1.0, 1.0/gamma )*refd
				else:
					dens  = refd
				# Append the area element
				parts.append([point,face[4],area,dens,cs])
		# Correct the area to fit to the real wall area
		fact = A / a
		for i in range(0,nu):
			for j in range(0,nv):
				parts[i*nv + j - nu*nv][2] = fact*parts[i*nv + j - nu*nv][2]
		idf   = idf + 1
		done  = min(idf*100 / len(faces), 100)
		if done != percentage:
			percentage = done
			stdout.write("%d%%..." % (percentage))
			stdout.flush()
	stdout.write("\n")
	return parts
