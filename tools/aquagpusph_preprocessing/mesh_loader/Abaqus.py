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
    end   = len(lines)
    for i in range(0,len(lines)):
        if(lines[i].upper().find('NODE') != -1):
            start = i+1
            break
    for i in range(i,len(lines)):
        # Remove all the spaces and tabulators
        line = lines[i].replace("\t", "").replace(" ", "").replace("\r", "").replace("\n", "")
        # Test if is not a point
        if (line == "") or (line.upper().find('ELEMENT') != -1):
            end = i
            break
    nPoints = end - start
    if nPoints <= 0:
        print("\tError: No points found")
        return None
    print("\t\tFound {0} points".format(nPoints))
    points  = {}
    i = start
    for i in range(start, end):
        # Remove all the spaces and tabulators
        line = lines[i].replace("\t", "").replace(" ", "").replace("\r", "").replace("\n", "")
        # Fields separation
        line = line.split(",")
        # Setup the point
        point = empty( (3) )
        for j in range(0,3):
            point[j] = float(line[j+1])
        points[line[0]] = point
    # Read the elements
    print("\tReading quadrilateral elements...")
    start = -1
    end   = len(lines)
    for i in range(0,len(lines)):
        if(lines[i].upper().find('ELEMENT') != -1) and (lines[i].upper().find('TYPE=S4R') != -1):
            start = i+1
            break
    for i in range(i,len(lines)):
        # Remove all the spaces and tabulators
        line = lines[i].replace("\t", "").replace(" ", "").replace("\r", "").replace("\n", "")
        # Test if is not a point
        if (line == "") or (line.upper().find('NODE') != -1):
            end = i
            break
    nFaces  = end - start
    if nFaces <= 0:
        print("\tError: No elements found")
        return None
    faces   = empty( (nFaces,5,3) )
    print("\t\tFound {0} faces".format(nFaces))
    faceID = 0
    for i in range(start,end):
        # Remove all the spaces and tabulators
        line = lines[i].replace("\t", "").replace(" ", "").replace("\r", "").replace("\n", "")
        # Fields separation
        line = line.split(",")
        # Read the points
        for j in range(0,4):
            faces[faceID,j,:] = points[line[j+1]][:]
        # Compute the normal
        vec1 = faces[faceID,2] - faces[faceID,0]
        vec2 = faces[faceID,3] - faces[faceID,1]
        n    = cross(vec1, vec2)
        faces[faceID,4] = n / sqrt( n[0]**2 + n[1]**2 + n[2]**2 )
        faceID = faceID + 1
    return faces
