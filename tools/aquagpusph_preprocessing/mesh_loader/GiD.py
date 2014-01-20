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

def loadQuadMesh(filename):
	""" Loads a GiD mesh file with quadrangular panels.
	@param filename Mesh file name.
	@return Faces defined by their 4 vertexes and a normal.
	"""
	f       = open(filename, 'r')
	lines   = f.readlines()
	# Read the points
	print("\tReading points...")
	start = -1
	end   = -1
	for i in range(0,len(lines)):
		if(lines[i].upper().find('COORDINATES') != -1):
			start = i+1
			break
	for i in range(i,len(lines)):
		if(lines[i].upper().find('END COORDINATES') != -1):
			end = i
	if (start < 0) or (end < start):
		print("\tError: No points found")
		return None
	nPoints = end - start
	points  = empty( (nPoints,3) )
	print("\t\tFound {0} points".format(nPoints))
	for i in range(0,nPoints):
		line  = lines[i+start].split()
		index = int(line[0])
		if(index != i+1):
			print("\t\tWarning: Point unsorted, {0} -> {1}".format(i+1, index))
		for j in range(0,3):
			points[index-1,j] = float(line[j+1])
	# Read the elements
	print("\tReading elements...")
	start = -1
	end   = -1
	for i in range(0,len(lines)):
		if(lines[i].upper().find('ELEMENTS') != -1):
			start = i+1
			break
	for i in range(i,len(lines)):
		if(lines[i].upper().find('END ELEMENTS') != -1):
			end = i
	if (start < 0) or (end < start):
		print("\tError: No elements found")
		return None
	nFaces  = end - start
	faces   = empty( (nFaces,5,3) )
	print("\t\tFound {0} faces".format(nFaces))
	for i in range(0,nFaces):
		line = lines[i+start].split()
		index = int(line[0])
		if(index != i+1):
			print("\t\tWarning: Elements unsorted, {0} -> {1}".format(i+1, index))
		if len(line) != 5:
			print("\t\tError: Element {0} is not quadrilateral".format(index))
			return None
		for j in range(0,4):
			faces[index-1,j,:] = points[int(line[j+1])-1,:]
		vec1 = faces[index-1,2] - faces[index-1,0]
		vec2 = faces[index-1,3] - faces[index-1,1]
		n    = cross(vec1, vec2)
		faces[index-1,4] = n / sqrt( n[0]**2 + n[1]**2 + n[2]**2 )
	return faces
