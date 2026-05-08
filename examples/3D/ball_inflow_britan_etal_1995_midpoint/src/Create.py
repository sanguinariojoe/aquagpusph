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
import platform
import math

# Constants & conditions
# ========================

# Chamber dimensions and ball radius
L = 0.6
B = 0.08
H = 0.08
R = 0.019

g = 0.0
hfac = 2.0
dr = 0.004

alpha = 0.0
delta = 1.0
visc_dyn = 0.0

courant = 0.25

gamma = 1.4

p1 = 1.0e5
rho1 = 1.0
c1 = math.sqrt(gamma * p1 / rho1)
M1 = 1.5
v1 = c1 * M1
u1 = 0.0

p2 = 2.45 * p1
rho2 = 1.86 * rho1
c2 = math.sqrt(gamma * p2 / rho2)
M2 = 0.70
v2 = c2 * M2
u2 = v1 - v2
print(f"Shock wave advancing at {v1} m/s")
print(f"Motion after shock wave {u2} m/s")

print(f"Sound speed before shock wave {c1} m/s")
print(f"Sound speed after shock wave {c2} m/s")

e1 = p1 / ((gamma - 1.0) * rho1)
e2 = p2 / ((gamma - 1.0) * rho2)

cs = max(c1, c2)

t_max = 1.5e-3

rhop = 89.4
refd = 0.0

sep = 2.0
h = hfac * dr

# Create the ball at the (0, 0, 0)
# ================================

# Get the number of subdivions of the icosphere so we get so far the same
# distance than dr
s = int(np.ceil(np.log(dr / R) / np.log(0.5)))
mesh = trimesh.creation.icosphere(subdivisions=s, radius=R)
print("Writing ball...")
n_ball = len(mesh.faces)
output = open("ball.dat", "w")
output.write("# r.x, r.y, r.z, r.w")
output.write(", normal.x, normal.y, normal.z, normal.w")
output.write(", tangent.x, tangent.y, tangent.z, tangent.w")
output.write(", u.x, u.y, u.z, u.w")
output.write(", dudt.x, dudt.y, dudt.z, dudt.w")
output.write(", rho, drhodt, e, dedt, m, imove\n")
for face, s, n in zip(mesh.faces, mesh.area_faces, mesh.face_normals):
    verts = mesh.vertices[face]
    r = np.mean(verts, axis=0)
    t = verts[1] - verts[0]
    t /= np.linalg.norm(t)
    n = -n
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
output.close()

# Create the fluid into the shock chamber
# =======================================

print("Writing fluid...")

Nx = int(round(L / dr))
Ny = int(round(B / dr))
Nz = int(round(H / dr))

L = Nx * dr
B = Ny * dr
H = Nz * dr

points = []
hL = 0.5 * L
hB = 0.5 * B
hh = 0.5 * H

x = np.linspace(-hL + 0.5 * dr, hL - 0.5 * dr, num=Nx)
y = np.linspace(-hB + 0.5 * dr, hB - 0.5 * dr, num=Ny)
z = np.linspace(-hh + 0.5 * dr, hh - 0.5 * dr, num=Nz)

xv, yv, zv = np.meshgrid(x, y, z)
points = np.asarray([xv.flatten(), yv.flatten(), zv.flatten()]).transpose()
print(f"{len(points)} candidate points")

# Exclude the points inside the sphere
distance = np.linalg.norm(points, axis=-1)
points = points[distance > R + 0.25 * dr]
print(f"{len(points)} points")

n_fluid = len(points)
output = open("fluid.dat", "w")
output.write("# r.x, r.y, r.z, r.w")
output.write(", normal.x, normal.y, normal.z, normal.w")
output.write(", tangent.x, tangent.y, tangent.z, tangent.w")
output.write(", u.x, u.y, u.z, u.w")
output.write(", dudt.x, dudt.y, dudt.z, dudt.w")
output.write(", rho, drhodt, e, dedt, m, imove\n")
for point in points:
    x, y, z = point
    imove = 1

    # The shock wave is placed just before the ball can sense it
    if x < -R - 0.5 * dr - sep * h:
        rho, ener, imove, velx = rho2, e2, 1, u2  
    else:
        rho, ener, imove, velx = rho1, e1, 1, u1

    mass = rho * dr**3.0
    string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}, {}, {}\n").format(
        x, y, z,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        velx, 0.0, 0.0,
        0.0, 0.0, 0.0,
        rho,
        0.0,
        ener,
        0.0,
        mass,
        imove)
    output.write(string)

# Create the shock chamber walls
# ==============================
# Bottom
for i in range(Nx):
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
            
            
# Right
for j in range(Ny):
    y = -hB + 0.5 * dr + j * dr
    for k in range(Nz):
        z = -hh + 0.5 * dr + k * dr
        x = hL
        nx = 1
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

#Setup the Inlet buffer particles. In this case we need to continuously feed
# with particles at a rate of Ny * Nz particles each dr / u2 seconds, during
# the full simulation. That is because we have no outlet, so the buffer will
# not be refilled during the runtime
ddom = 2 * sep * h
domain_min = (-hL - ddom, -hB - ddom, -hh - ddom, 0.0)
domain_max = (hL + ddom, hB + ddom, hh + ddom, 0.0)

n_buffer_depth = 16
n_buffer = Ny * Nz * (n_buffer_depth + int(math.ceil(u2 / dr * t_max)))
print(f"{n_buffer} buffer particles")

x = domain_max[0] + sep * h
y = domain_max[1] + sep * h
z = domain_max[2] + sep * h
for i in range(n_buffer):
    #n += 1
    imove = -255       
    mass = rho2 * dr**2.0
    string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}, {}, {}\n").format(
        x, y, z,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
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

domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(rho1),
        'RHO1':str(rho1), 'RHO2':str(rho2), 'RHOP':str(rhop),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'L':str(L), 'B':str(B), 'H':str(H), 'R':str(R), 'GAMMA':str(gamma),        
        'NX':str(Nx), 'NY':str(Ny), 'NZ':str(Nz), 'T': str(t_max),
        'U2':str(u2), 'E2':str(e2),'L_2':str(L/2.0), 'B_2':str(B/2.0), 'H_2':str(H/2.0), 
        'n_ball':str(n_ball), 'n_fluid':str(n_fluid)}
exttool_lib_name = "britan_ball_sim_com_inflow.dll" if platform.system() == "Windows" \
    else "libbritan_ball_sim_com_inflow.so"
data['EXTTOOL_LIB_PATH'] = os.path.join(script_folder, exttool_lib_name)
utils.configure(data, os.path.join(script_folder, "templates"))
