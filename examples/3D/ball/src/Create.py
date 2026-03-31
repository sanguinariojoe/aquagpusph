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

import os
import sys
script_folder = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_folder, "../../"))
import aqua_example_utils as utils

import numpy as np
import trimesh
import meshio
import platform
import math

# Constants & conditions
# ========================

g = 0.0
hfac = 2.0
dr = 0.004

#refd = 1.0
alpha = 0.0
delta = 1.0
visc_dyn = 0.0

courant = 0.25

gamma = 1.4

p1 = 1.0e5
rho1 = 1.0
c1 = math.sqrt(gamma * p1 / rho1)
M1=1.5
v1=c1*M1
u1=0.0

p2=2.45*p1
rho2=1.86*rho1
c2 = math.sqrt(gamma * p2 / rho2)
M2=0.70
v2=c2*M2
u2=v1-v2

e1 = p1 / ((gamma - 1.0) * rho1)
e2 = p2 / ((gamma - 1.0) * rho2)

cs = max(c1, c2)

t_max=1.5e-3

rhop=89.4

L = 0.6
B = 0.08
H = 0.08
R = 0.019

sep = 2.0
h = hfac * dr
Lext = L

TO_METERS = 1.e-3
IN_POINT = [0, 0, 0]

# Read the surface mesh file and create the boundary particles
# ============================================================
files = [f for f in os.listdir(script_folder) if os.path.isfile(
    os.path.join(script_folder, f))]
files = sorted(
    [f for f in files if f.startswith('ball.') and f.endswith('.stl')])
files = sorted(
    [f for f in files if not f.endswith('subdivided.stl')])
n_balls = len(files)

meshes = []
xml_files = []
prefixes = []
bbox = np.array([[np.finfo(np.float64).max,
                  np.finfo(np.float64).max,
                  np.finfo(np.float64).max],
                 [np.finfo(np.float64).min,
                  np.finfo(np.float64).min,
                  np.finfo(np.float64).min]])

for i, f in enumerate(files):
    fout = f[:-4] + ".subdivided.stl"
    mesh = trimesh.load(os.path.join(script_folder, f))
    meshes.append(mesh)
    new_v, new_f = trimesh.remesh.subdivide_to_size(mesh.vertices,
                                                    mesh.faces,
                                                    dr)
    mesh = trimesh.Trimesh(vertices=new_v, faces=new_f)
    mesh.export(fout)
    mesh = meshio.read(os.path.join(script_folder, fout))
    fout = f[:-4] + ".dat"
    print(f"Writing {fout}...")
    output = open(f"{fout}", "w")
    output.write("# r.x, r.y, r.z, r.w")
    output.write(", normal.x, normal.y, normal.z, normal.w")
    output.write(", tangent.x, tangent.y, tangent.z, tangent.w")
    output.write(", u.x, u.y, u.z, u.w")
    output.write(", dudt.x, dudt.y, dudt.z, dudt.w")
    output.write(", rho, drhodt, e, dedt, m, imove\n")
    n_parts = 0
    verts = mesh.points
    for cell in mesh.cells:
        def triangle(cell, verts):
            a, b, c = [verts[i] for i in cell]
            r = np.mean([a, b, c], axis=0)
            t = b - a
            n = np.cross(b - a, c - a)
            s = 0.5 * np.linalg.norm(n)
            t /= np.linalg.norm(t)
            n /= 2.0 * s
            n = -n  # Inwards normals
            return r, n, t, s

        for elem in cell.data:
            r, n, t, s = triangle(elem, verts)
            bbox[0] = np.min((bbox[0], r), axis=0)
            bbox[1] = np.max((bbox[1], r), axis=0)
            dens = refd
            imove = -3
            string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}, {}, {}\n").format(
                r[0], r[1], r[2],
                n[0], n[1], n[2],
                t[0], t[1], t[2],
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                rho2,
                0.0,
                e2,
                0.0,
                s,
                imove)

            output.write(string)
            n_parts += 1
    output.close()

    prefix = fout[:-4].replace('.', '_') + '_'
    prefixes.append(prefix)
    data = {'N_PARTS':str(n_parts), 'REFD':str(rhop),
            'VISC_DYN':str(visc_dyn), 'DELTA':str(delta),
            'FIN':fout, 'FOUT':fout[:-4], 'PREFIX':prefix, 'ISET':str(i + 1)}
    utils.configure(data, os.path.join(script_folder, "ball_template"))
    fout = f[:-4] + ".xml"
    print(f"Writing {fout}...")
    os.rename('ball.xml', fout)
    xml_files.append(fout)

# Write a general reader to be included from Main.xml
with open("ball.xml", "w") as f:
    f.write("<sphInput>\n")
    for xml in xml_files:
        f.write(f'\t<Include file="{xml}" />\n')
    f.write("\t<Tools>\n")
    for prefix in prefixes[1:]:
        f.write(f'\t\t<Tool name="{prefix}cfd BIe backup p" action="try_remove" type="dummy"></Tool>\n')
        f.write(f'\t\t<Tool name="{prefix}cfd BIe backup force_visc" action="try_remove" type="dummy"></Tool>\n')
        f.write(f'\t\t<Tool name="{prefix}cfd BIe backup moment_visc" action="try_remove" type="dummy"></Tool>\n')
    f.write("\t</Tools>\n")
    f.write("</sphInput>\n")

# Fluid
# ============

print("Writing fluid...")

Nx = nx = int(round(Lext / dr))
Ny = ny = int(round(B / dr))
Nz = nz = int(round(H / dr))


Lext = Nx * dr
L = Nx * dr
B = Ny * dr
H = Nz * dr

points = []
hL = 0.5 * Lext
hB = 0.5 * B
hh = 0.5 * H

x = np.linspace(-hL + 0.5 * dr, hL - 0.5 * dr, num=Nx)
y = np.linspace(-hB + 0.5 * dr, hB - 0.5 * dr, num=Ny)
z = np.linspace(-hh + 0.5 * dr, hh - 0.5 * dr, num=Nz)

xv, yv, zv = np.meshgrid(x, y, z)
points = np.asarray([xv.flatten(), yv.flatten(), zv.flatten()]).transpose()
print(f"{len(points)} candidate points")


# Removing points in the sphere
# ================================
for i, mesh in enumerate(meshes):
    [[xmin, ymin, zmin], [xmax, ymax, zmax]] = mesh.bounds
    mask = (points[:, 0] >= (xmin - 0.5 * dr)) & \
           (points[:, 0] <= (xmax + 0.5 * dr)) & \
           (points[:, 1] >= (ymin - 0.5 * dr)) & \
           (points[:, 1] <= (ymax + 0.5 * dr)) & \
           (points[:, 2] >= (zmin - 0.5 * dr)) & \
           (points[:, 2] <= (zmax + 0.5 * dr))
    mask[mask] = np.asarray(mesh.contains(points[mask]))
    print(f"Dropping {np.sum(mask)} points inside {xml_files[i][:-4]}")
    points = points[np.logical_not(mask)]
    distance = dr * np.ones(len(points))
    mask = (points[:, 0] >= (xmin - 0.5 * dr)) & \
           (points[:, 0] <= (xmax + 0.5 * dr)) & \
           (points[:, 1] >= (ymin - 0.5 * dr)) & \
           (points[:, 1] <= (ymax + 0.5 * dr)) & \
           (points[:, 2] >= (zmin - 0.5 * dr)) & \
           (points[:, 2] <= (zmax + 0.5 * dr))
    _, distance[mask], _ = mesh.nearest.on_surface(points[mask])
    mask = distance < 0.25 * dr
    print(f"Dropping {np.sum(mask)} points too close to {xml_files[i][:-4]}")
    points = points[np.logical_not(mask)]

# Go for the remaining points
# =============================

output = open("fluid.dat", "w")
output.write("# r.x, r.y, r.z, r.w")
output.write(", normal.x, normal.y, normal.z, normal.w")
output.write(", tangent.x, tangent.y, tangent.z, tangent.w")
output.write(", u.x, u.y, u.z, u.w")
output.write(", dudt.x, dudt.y, dudt.z, dudt.w")
output.write(", rho, drhodt, e, dedt, m, imove\n")
n_fluid = 0

for point in points:
    x, y, z = point
    imove = 1
    
    mod=np.sqrt(x*x+y*y+z*z)
    
    if x < -R -0.5*dr:
        rho, ener, imove = rho2, e2, 1  
        
    else:
        rho, ener, imove = rho1, e1, 1
        
    mass = rho * dr**3.0
    string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}, {}, {}\n").format(
        x, y, z,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        rho,
        0.0,
        ener,
        0.0,
        mass,
        imove)
    output.write(string)
    n_fluid += 1

# Bottom
for i in range(Nx):
    #x = -hL - 0.5 * dr + i * dr
    x = -hL + 0.5 * dr + i * dr
    for j in range(Ny):
        y = -hB + 0.5 * dr + j * dr
        z = -hh
        imove = -3
        mass = dr**2.0
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}, {}, {}\n").format(
            x, y, z,
            0.0, 0.0, -1.0,
            -1.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            rho2,
            0.0,
            e2,
            0.0,
            mass,
            imove)
        output.write(string)
        n_fluid += 1

# Top
for i in range(Nx):
    #x = -hL - 0.5 * dr + i * dr
    x = -hL + 0.5 * dr + i * dr
    for j in range(Ny):
        y = -hB + 0.5 * dr + j * dr
        z = hh
        imove = -3
        mass = dr**2.0
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}, {}, {}\n").format(
            x, y, z,
            0.0, 0.0, 1.0,
            1.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            rho2,
            0.0,
            e2,
            0.0,
            mass,
            imove)
        output.write(string)
        n_fluid += 1

# Front and back
for i in range(Nx):
    #x = -hL - 0.5 * dr + i * dr
    x = -hL + 0.5 * dr + i * dr
    for k in range(Nz):
        z = -hh + 0.5 * dr + k * dr
        for j in (-1, 1):
            y = hB * j
            ny = j
            imove = -3
            mass = dr**2.0
            string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}, {}, {}\n").format(
                x, y, z,
                0.0, ny, 0.0,
                0.0, 0.0, -ny,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                rho2,
                0.0,
                e2,
                0.0,
                mass,
                imove)
            output.write(string)
            n_fluid += 1
            
            
# two sides        
for j in range(Ny):
    y = -hB + 0.5 * dr + j * dr
    for k in range(Nz):
        z = -hh + 0.5 * dr + k * dr
        for i in (-1, 1):
            x = hL * i #+ dr * i
            nx = i
            imove = -3
            mass = dr**2.0
            string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}, {}, {}\n").format(
                x, y, z,
                nx, 0.0, 0.0,
                0.0, 0.0, -nx,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                rho2,
                0.0,
                e2,
                0.0,
                mass,
                imove)
            output.write(string)
            n_fluid += 1
            
            
output.close()

domain_min = (-Lext, -B, -0.5 * H, 0.0)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (Lext, B, 1.5 * H, 0.0)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(rho1),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'L':str(L), 'B':str(B), 'H':str(H), 'GAMMA':str(gamma),        
        'NX':str(Nx), 'NY':str(Ny), 'NZ':str(Nz), 'T': str(t_max),
        'NROCKS':str(1), 'n_fluid':str(n_fluid)}
exttool_lib_name = "rocks_sim.dll" if platform.system() == "Windows" \
    else "libball_sim.so"
data['EXTTOOL_LIB_PATH'] = os.path.join(script_folder, exttool_lib_name)
utils.configure(data, os.path.join(script_folder, "templates"))
