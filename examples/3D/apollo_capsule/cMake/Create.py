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

import os.path as path
import numpy as np
from scipy.spatial.transform import Rotation
import meshio

# Input data
# ==========

COG = [0.142, 0.0, 0.318]
PITCH = -12.0
VEL = 9.88
DT = 1e-6
Ma = 1e-1  # This is just used for the variable time step
courant = 0.1
g = 9.81
hfac = 2.0
dr = 0.16
cs = 1000.0
refd = 998.0
alpha = 0.0
visc_dyn = 0.000894
# Tank dimensions
H = 16.0
L = 16.0
# Fluid
h = 8.0

# Estimate the simulation time
z0 = 2.0 * hfac * dr
assert VEL**2 - 2 * g * z0 > 0.0
v0 = -np.sqrt(VEL**2 - 2 * g * z0)
t0 = (v0 + VEL) / g
VEL = -v0
T = t0 + 0.18
print(f"Releasing velocity = {VEL} m / s")
print(f"Flying time = {t0} s")

# Read the mesh file and create the apollo capsule particles
# =========================================================
mesh = meshio.read("CAD/Mesh.med")
verts = mesh.points
Rot = Rotation.from_euler('y', PITCH, degrees=True)
verts = Rot.apply(verts)
cog = Rot.apply(COG)
verts = np.subtract(verts, cog)
minz = np.min(verts[:, -1])
cogz = z0 - minz
verts[:, -1] += cogz
maxz = np.max(verts[:, -1])
print(f"Initial COGz = {cogz} m")
print("Writing apollo...")
output = open("Titan.dat", "w")
output.write("# r.x, r.y, r.z, r.w")
output.write(", normal.x, normal.y, normal.z, normal.w")
output.write(", tangent.x, tangent.y, tangent.z, tangent.w")
output.write(", u.x, u.y, u.z, u.w")
output.write(", dudt.x, dudt.y, dudt.z, dudt.w")
output.write(", rho, drhodt, m, imove\n")
n_apollo = 0
for cell in mesh.cells:
    def triangle(cell, verts):
        a, b, c = [verts[i] for i in cell]
        r = np.mean([a, b, c], axis=0)
        t = b - a
        n = np.cross(b - a, c - a)
        s = 0.5 * np.linalg.norm(n)
        t /= np.linalg.norm(t)
        n /= 2.0 * s
        if np.dot(r - [0, 0, cogz], n) > 0.0:
            n = -n
        return r, n, t, s

    def quad(cell, verts):
        a, b, c, d = [verts[i] for i in cell]
        r = np.mean([a, b, c, d], axis=0)
        t = b - a
        t /= np.linalg.norm(t)
        n = np.cross(c - a, d - b)
        n /= np.linalg.norm(n)
        if np.dot(r - [0, 0, cogz], n) > 0.0:
            n = -n
        s1 = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
        s2 = 0.5 * np.linalg.norm(np.cross(b - d, c - d))
        s = s1 + s2
        return r, n, t, s
        
    if cell.type not in ['quad', 'triangle']:
        continue
    backend = quad if cell.type == 'quad' else triangle
    for elem in cell.data:
        r, n, t, s = backend(elem, verts)
        dens = refd
        imove = -3
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            r[0], r[1], r[2],
            n[0], n[1], n[2],
            t[0], t[1], t[2],
            0.0, 0.0, -VEL,
            0.0, 0.0, -g,
            dens,
            0.0,
            s,
            imove)
        output.write(string)
        n_apollo += 1
output.close()

print("Writing pool...")
Ny = ny = Nx = nx = int(round(L / dr))
Nz = int(round(H / dr))
nz = int(round(h / dr))
output = open("Pool.dat", "w")
output.write("# r.x, r.y, r.z, r.w")
output.write(", normal.x, normal.y, normal.z, normal.w")
output.write(", tangent.x, tangent.y, tangent.z, tangent.w")
output.write(", u.x, u.y, u.z, u.w")
output.write(", dudt.x, dudt.y, dudt.z, dudt.w")
output.write(", rho, drhodt, m, imove\n")
n_pool = 0
for i in range(nx):
    x = -0.5 * L + 0.5 * dr + i * dr
    for j in range(ny):
        y = -0.5 * L + 0.5 * dr + j * dr
        for k in range(nz):
            z = -h + 0.5 * dr + k * dr
            imove = 1
            press = -refd * g * z
            dens = refd + press / cs**2 
            mass = dens * dr**3.0
            string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
                x, y, z,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                dens,
                0.0,
                mass,
                imove)
            output.write(string)
            n_pool += 1

# Bottom
for i in range(Nx):
    x = -0.5 * L + 0.5 * dr + i * dr
    for j in range(Ny):
        y = -0.5 * L + 0.5 * dr + j * dr
        z = -h
        imove = -3
        press = -refd * g * z
        dens = refd + press / cs**2 
        mass = dr**2.0
        string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
            x, y, z,
            0.0, 0.0, -1.0,
            -1.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            dens,
            0.0,
            mass,
            imove)
        output.write(string)
        n_pool += 1

# Front and back
for i in range(Nx):
    x = -0.5 * L + 0.5 * dr + i * dr
    for k in range(Nz):
        z = -h + 0.5 * dr + k * dr
        for j in (-1, 1):
            y = 0.5 * L * j
            ny = j
            imove = -3
            press = -refd * g * z
            dens = refd + press / cs**2 
            mass = dr**2.0
            string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
                x, y, z,
                0.0, ny, 0.0,
                0.0, 0.0, -ny,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                dens,
                0.0,
                mass,
                imove)
            output.write(string)
            n_pool += 1

# Left and right
for j in range(Ny):
    y = -0.5 * L + 0.5 * dr + j * dr
    for k in range(Nz):
        z = -h + 0.5 * dr + k * dr
        for i in (-1, 1):
            x = 0.5 * L * i
            nx = i
            imove = -3
            press = -refd * g * z
            dens = refd + press / cs**2 
            mass = dr**2.0
            string = ("{} {} {} 0.0, " * 5 + "{}, {}, {}, {}\n").format(
                x, y, z,
                nx, 0.0, 0.0,
                0.0, 0.0, -ny,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                dens,
                0.0,
                mass,
                imove)
            output.write(string)
            n_pool += 1
output.close()

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('Fluids.xml', 'Main.xml', 'Motion.xml', 'Settings.xml', 'SPH.xml',
       'Time.xml', 'plot_m.py', 'plot_e.py')

domain_min = (-L, -L, -H, 0.0)
domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = (L, L, float(maxz + 4.0 * hfac * dr), 0.0)
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'G':str(g), 'MA':str(Ma), 'DTMIN':str(DT),
        'COGZ':str(cogz), 'VEL':str(VEL), 'PITCH':str(PITCH),
        'T':str(T), 'n_apollo':str(n_apollo), 'n_pool':str(n_pool)}
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
