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
from .vec import *
import sys
from sys import stdout
import os
import time
import multiprocessing
from multiprocessing import Process, Queue

def bounds(faces):
    """ Look for solid bounds.
    @param faces Faces defined by 4 points and a normal
    @return two points as the minimum and maximum coordinates.
    """
    minimum = copy(faces[0][0])
    maximum = copy(faces[0][0])
    for f in faces:
        for i in range(0,4):
            minimum[0] = min(minimum[0], f[i][0])
            minimum[1] = min(minimum[1], f[i][1])
            minimum[2] = min(minimum[2], f[i][2])
            maximum[0] = max(maximum[0], f[i][0])
            maximum[1] = max(maximum[1], f[i][1])
            maximum[2] = max(maximum[2], f[i][2])
    return [minimum, maximum]

def isPointInWall(p, w, threshold=0.001):
    """ Test if a point is into a wall.
    @param p Point to test.
    @param w Wall to test (defined by 3 vertices).
    @param threshold If the point is not 
    at triangle plane, but is nearby this
    value, can be considereed valid.
    @return True if the point is inside
    the wall, False otherwise.
    @note Use baricentric algorithm.
    """
    V0 = w[1] - w[0]
    V1 = w[2] - w[0]
    n  = cross(V0,V1)
    N  = linalg.norm(n)
    if not N:
        return False
    n  = n/N
    V2 = p - w[0]
    dist = dot(V2,n)
    if abs(dist) >= threshold:    # Is not in the plane
        return False
    V2 = V2 - n*dist
    # Compute dots
    dot00 = dot(V0, V0)
    dot01 = dot(V0, V1)
    dot02 = dot(V0, V2)
    dot11 = dot(V1, V1)
    dot12 = dot(V1, V2)
    # Get triangle UV coordinates
    div = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * div
    v = (dot00 * dot12 - dot01 * dot02) * div
    # Check valid points
    if u >= 0.0 and v >= 0.0 and u+v <= 1.0:
        return True
    return False

def isPointInSolid(faces, p, threshold=0.001):
    """ Test if a point is into a solid face.
    @param faces Faces defined by 4 points and a normal
    @param p Point.
    @return True if the point is on the solid surface
    False otherwise.
    """
    for f in faces:
        # Divide the face in two triangles
        walls = empty( (2,3,3) )
        walls[0][0] = f[0]
        walls[0][1] = f[1]
        walls[0][2] = f[2]
        walls[1][0] = f[0]
        walls[1][1] = f[2]
        walls[1][2] = f[3]
        if(isPointInWall(p,walls[0],threshold)):
            return True
        if(isPointInWall(p,walls[1],threshold)):
            return True
    return False

def getWallDist(p, w, maxDist=None):
    """ Return the distance to a wall, and the normal oriented to the point.
    @param p Point to test.
    @param w Wall to test (defined by 3 vertices).
    @return (flag,d,n)
    flag = True if the particle projected over the wall is inside it, False otherwise
    d    = Distance to the wall
    n    = Normal to the wall oriented to the point
    @note Use baricentric algorithm.
    """
    # Compute missoriented normal
    V0 = w[1] - w[0]
    V1 = w[2] - w[0]
    n = cross(V0,V1)
    N = linalg.norm(n)
    if not N:
        return (False,0,zeros(3))
    n = n/N
    # Compute the distance to the wall
    V2 = p - w[0]
    d = dot(V2,n)
    # Swap the normal if needed
    if d < 0.0:
        d = -d
        n = -n
    # Test if we are not near enough
    if maxDist is not None and maxDist < d:
        return False,d,n
    # Test if the particle projected over the wall is inside it
    flag = False
    V2 = V2 - n*d
    # Compute dots
    dot00 = dot(V0, V0)
    dot01 = dot(V0, V1)
    dot02 = dot(V0, V2)
    dot11 = dot(V1, V1)
    dot12 = dot(V1, V2)
    # Get triangle UV coordinates
    div = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * div
    v = (dot00 * dot12 - dot01 * dot02) * div
    # Check valid points
    if u >= 0.0 and v >= 0.0 and u+v <= 1.0:
        flag = True
    return flag,d,n

def fixSolidEffect(faces, p, destroyDist=0.001, fix=[]):
    """ Test if a point is too near to a solid wall (destroying it), or
    if is affected by the solid and may be moved to be to a correct distance
    to the wall.
    @param faces Faces defined by 4 points and a normal
    @param p Point.
    @param destroyDist Distance to destroy the particle (too near).
    @param fix Intervals of fixing (each inteval must be a tuple (d0,d1,d)).
    d0 = minimum distance.
    d1 = maximum distance.
    d  = selected distance.
    @return A suggested point if the particle can be better placed, or
    None if must be discarded.
    """
    p = copy(p)
    for f in faces:
        # Divide the face in two triangles
        walls = empty( (2,3,3) )
        walls[0][0] = f[0]
        walls[0][1] = f[1]
        walls[0][2] = f[2]
        walls[1][0] = f[0]
        walls[1][1] = f[2]
        walls[1][2] = f[3]
        for w in walls:
            flag,d,n = getWallDist(p, w, destroyDist)
            if not flag:
                continue
            if d < destroyDist:
                return None
            for f in fix:
                if (d < f[1]) and (d > f[0]):
                    p = p + n*(f[2] - d)
    return p

def intersection(p, v, w):
    """ Computes intersection of a ray with triangular wall.
    @param p Ray launching point.
    @param v Ray direction vector (must be normalized).
    @param w Triangular wall vertices.
    @return Intersection point, None if not exist.
    """
    # wall normal (not relevant if inner or outer)
    V0 = w[1] - w[0]
    V1 = w[2] - w[0]
    V2 = w[2] - w[1]
    n  = cross(V0,V1)
    N  = length(n)
    if not N:
        return None
    n  = n/N
    # Ensure that plane and line can be intersected
    if abs(dot(v,n)) < 0.001:
        return None
    # Get line-plane intersection
    d = dot(w[0] - p, n) / dot(v,n)
    point = p + v*d
    # Ensure that the point is into the triangle
    if not isPointInWall(point,w):
        return None
    return point

def intersections(faces, p, v):
    """ Computes all intersections of a ray with solid boundaries.
    @param faces Faces defined by 4 points and a normal
    @param p Ray launching point.
    @param v Ray vector.
    @return Array of intersection points.
    """
    points = []
    for f in faces:
        # Divide the face in two triangles
        walls = empty( (2,3,3) )
        walls[0][0] = f[0]
        walls[0][1] = f[1]
        walls[0][2] = f[2]
        walls[1][0] = f[0]
        walls[1][1] = f[2]
        walls[1][2] = f[3]
        point = intersection(p,v,walls[0])
        if point != None:
            points.append(point)
        point = intersection(p,v,walls[1])
        if point != None:
            points.append(point)
    return points

def filterFacesZ(z,dz,faces):
    """ Compute the list of faces that must be considered by a
    z coordinate, meaning that are present on this coordinate.
    @param z z coordinate
    @param dz Differential of height
    @param faces Faces defined by 4 points and a normal
    """
    filtered = []
    for f in faces:
        minZ = f[0][2]
        maxZ = f[0][2]
        for i in range(1,4):
            minZ = min(minZ, f[i][2])
            maxZ = max(maxZ, f[i][2])
        minZ = minZ - dz
        maxZ = maxZ + dz
        if(minZ < z and maxZ > z):
            filtered.append(f)
    return filtered

def filterFacesY(y,dy,faces):
    """ Compute the list of faces that must be considered by a
    y coordinate, meaning that are present on this coordinate.
    @param y y coordinate
    @param dy Differential of length
    @param faces Faces defined by 4 points and a normal
    """
    filtered = []
    for f in faces:
        minY = f[0][1]
        maxY = f[0][1]
        for i in range(1,4):
            minY = min(minY, f[i][1])
            maxY = max(maxY, f[i][1])
        minY = minY - dy
        maxY = maxY + dy
        if(minY < y and maxY > y):
            filtered.append(f)
    return filtered

def filterFacesX(x,dx,faces):
    """ Compute the list of faces that must be considered by a
    x coordinate, meaning that are present on this coordinate.
    @param x x coordinate
    @param dx Differential of length
    @param faces Faces defined by 4 points and a normal
    """
    filtered = []
    for f in faces:
        minX = f[0][1]
        maxX = f[0][1]
        for i in range(1,4):
            minX = min(minX, f[i][1])
            maxX = max(maxX, f[i][1])
        minX = minX - dx
        maxX = maxX + dx
        if(minX < x and maxX > x):
            filtered.append(f)
    return filtered

def zThread(q, z, dz, bbox, faces, level, dr, refd, cs, g):
    """ Creates the particles for a desired z level
    @param q Queue to use to set the particles
    @param z z coordinate to fill with particles
    @param dz Differential of height
    @param bbox Solid boundary box
    @param faces Faces defined by 4 points and a normal
    @param level Fluid level
    @param dr Distance between particles
    @param refd Fluid density
    @param cs Sound speed
    @param g Gravity acceleration
    """
    parts = []
    # To reduce the computational requirements we are filtering the faces used
    # to compute this z level
    facesZ = filterFacesZ(z,dz,faces)
    # Now we are on the z coordinate that we want to fill with fluid
    # particles, so we will start running on y coordinate launching rays on
    # the x direction looking for valid intervals into the solid mesh.
    Dy = bbox[1][1] - bbox[0][1]
    ny = int(round(Dy / dr))
    dy = Dy / ny
    y  = bbox[0][1] + 0.5*dy
    for i in range(0,ny):
        # We can filter the faces for this y coordinate as well in order to
        # reduce the computational efforts
        facesY = filterFacesY(y,dy,facesZ)
        # Now we can launch a ray on x direction looking for the intersections
        # with the solid faces. An even number of intersections may be
        # expected for closed solids, then the odd intervals (which is a pair
        # of intersection points) are the ones which particles can be placed
        p = empty(3)
        p[0] = bbox[0][0] - dr
        p[1] = y
        p[2] = z
        v = zeros(3)
        v[0] = 1.0
        intervals = intersections(facesY,p,v)
        if not len(intervals):
            y = y + dy
            continue
        # Sort the points (and remove repeated)
        aux = []
        for interval in intervals:
            aux.append(interval[0])
        aux = set(aux)
        intervals = []
        for interval in aux:
            intervals.append(interval)
        intervals.sort()
        # Ensure that we get a even number of intersections
        if len(intervals) % 2:
            print("ERROR: Solid seems to be unclosed.")
            print("\tIntersection with {} points found (odd value is unadmissible) in z = {}".format(len(intervals), z))
            q.put([])
            return
        j = 1
        while(j < len(intervals)):
            Dx = intervals[j] - intervals[j-1]
            nx = int(round(Dx / dr))
            for k in range(0,nx):
                dx = Dx / nx
                x = intervals[j-1] + (0.5 + k)*dx
                point = empty(3)
                normal = zeros(3)
                point[0] = x
                point[1] = y
                point[2] = z
                vol = dx*dy*dz
                pdr = vol**(1.0/3.0)
                # Discard the particles excesively near to the solid, or move
                # it to be in 0.5 dr or 1.5 dr.
                point = fixSolidEffect(facesY, point, 0.1*pdr,
                                       [[0.1*pdr,1.0*pdr,0.5*pdr],
                                        [1.0*pdr,2.0*pdr,1.5*pdr]])
                # Discard the particles excesively near to the solid (0.5*dr).
                # point = fixSolidEffect(facesY, point, 0.5*pdr)
                if point == None:
                    continue
                press = refd * g * (level-z)
                dens = refd + press / cs**2
                mass = dens * vol
                # Append the new particle
                parts.append([point,normal,mass,dens,cs])
            j = j+2
        y = y + dy
    q.put(parts)

def perform(faces, level, dr, refd, cs, g):
    """ Create fluid particles
    @param faces Faces defined by 4 points and a normal
    @param level Fluid level
    @param dr Distance between particles
    @param refd Fluid density
    @param cs Sound speed
    @param g Gravity acceleration
    @return Fluid particles
    """
    parts  = []
    nCores = multiprocessing.cpu_count()
    # To perform particles we will start at minimum z coordinate of
    # the solid, and we will start going up until we get the level
    # required. Since the specified distance between particles can
    # cause that level will not be right, so must be corrected in
    # order to fit it.
    bbox = bounds(faces)
    if(level <= bbox[0][2]):
        print("WARNING: level is not enough to create any fluid particles")
        print("\tMinimum solid height is %g, but %g level has been provided" % (bbox[0][2], level))
        return parts
    if(level > bbox[1][2]):
        print("WARNING: level is greather than the maximum solid height")
        print("\tLevel will be corrected therefore to %g meters" % (bbox[1][2]))
        level = bbox[1][2]
    Dz   = level - bbox[0][2]
    nz   = int(round(Dz / dr))
    dz   = Dz / nz
    z    = bbox[0][2] + 0.5*dz
    idz  = 0
    # In order to process particles in parallel we will create as many threads
    # as available on system
    q           = []
    p           = []
    nValidCores = 0
    for i in range(0, nCores):
        q.append(None)
        p.append(None)
    percentage = 0
    stdout.write("{}%...".format(percentage))
    stdout.flush()
    # We will iterate until the threads inside the loop can end their work
    while(True):
        # Check if we have end the work
        if(z >= level and nValidCores == 0):
            # We are waiting for all remaining threads ends their work
            time.sleep(0.1)
            break
        # Check for unactive threads and assign them new work
        for i in range(0, nCores):
            if(z >= level):
                # No new work to assign
                break
            if( q[i] != None and p[i] != None ):
                continue
            q[i] = Queue()
            p[i] = Process(target=zThread, args=(q[i],z,dz,bbox,faces,level,dr,refd,cs,g))
            p[i].start()
            nValidCores = nValidCores + 1
            z = z + dz
        # Check if some thread is returning values
        for i in range(0, nCores):
            if( q[i] == None and p[i] == None ):
                continue
            if(not q[i].empty()):
                # Then a thread is returning particles, get it and destroy
                parts.extend(q[i].get())
                p[i].join()
                q[i]        = None
                p[i]        = None
                nValidCores = nValidCores - 1
                idz         = idz + 1
                percentage  = min(idz*100 / nz, 100)
                stdout.write("{}%...".format(percentage))
                stdout.flush()
        # Let a little bit of time to the threads in roder to don't saturate this one
        time.sleep(0.1)
    stdout.write("\n")
    stdout.flush()
    return parts
