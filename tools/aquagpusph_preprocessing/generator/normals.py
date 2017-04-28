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

import numpy as np
from .fluid import intersections, isPointInSolid


def is_odd(num):
   return num % 2 != 0


def recompute_normals(faces):
    """Recompute the normals to assert they are outward oriented.

    Position arguments:
    faces -- List of faces. The faces will be overwritten, the user is
    responsible to make a copy before calling this method in case the original
    normals should be preserved.

    Returned value:
    List of faces with the normals recomputed.
    """
    for i,face in enumerate(faces):
        p = np.sum(face[0:len(face) - 1], axis=0) / 4.0
        n = face[-1]
        # Get the intersections relative to the face center
        points = np.add(np.asarray(intersections(faces, p, n)), -p)
        # Get the votes
        dists = np.dot(points, n)
        try:
            noswap = len(np.unique(dists[np.where(dists > 0.001)]))
        except TypeError:
            noswap = 0
        try:
            swap = len(np.unique(dists[np.where(dists < -0.001)]))
        except TypeError:
            swap = 0
        # Check for errors
        if is_odd(noswap) and is_odd(swap):
            print("ERROR: Odd number of intersections found at both directions of the face {}".format(i))
            continue
        if (not is_odd(noswap)) and (not is_odd(swap)):
            print("ERROR: Even number of intersections found at both directions of the face {}".format(i))
            continue
        # Swap the normal if required
        if not is_odd(swap):
            print("Changing normal orientation of the face {}".format(i))
            faces[i][4] *= -1.0
    return faces
