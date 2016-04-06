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
import math

# Input data
# ==========

g = 0.0
hfac = 4.0
courant = 0.2
refd = 1.0
alpha = 0.0
delta = 1.0
U = 1.0
cs = 20.0 * U
Re = 100.0
p0 = 3.0 * refd * U**2
BC = 'GP'  # BI or GP
# Cylinder and Channel dimensions
D = 1.0
L = 10.0 * D
H = 5.0 * D
# Position of the cylinder
x_square = 0.15 * L
y_square = 0.5 * H
# Number of fluid particles in y direction
ny = 400

# Distance between particles
# ==========================
sep = 2.0
dr = H / ny
h = hfac * dr
domain_min = (-3.0 * sep * h, -3.0 * sep * h)
domain_max = (L + 3.0 * sep * h, H + 3.0 * sep * h)

# Artificial viscosity
# ====================
visc_dyn = refd * U * D / Re
visc_dyn = max(alpha / 8.0 * refd * hfac * dr * cs, visc_dyn)

# Particles generation
# ====================
def write_particle(fhandler, r, n=(0,0), u=(0,0), dudt=(0,0),
                   rho=refd, drhodt=0, mass=refd*dr**2, imove=1,
                   associated='N'):
    string = ("{} {}, " * 4 + "{}, {}, {}, {}, {}\n").format(
            r[0], r[1],
            n[0], n[1],
            u[0], u[1],
            dudt[0], dudt[1],
            rho,
            drhodt,
            mass,
            imove,
            associated)
    fhandler.write(string)


def write_gp_particles(fhandler, r, n, associated, m_fact=1.0):
    """Write the set of Ghost particles behind an area boundary element.
    
    Position arguments:
    r -- Position of the area element
    n -- Outward normal of the area element
    associated -- Index of the mirroring particle (the area element)
    
    Keyword arguments:
    m_fact -- Mass multiplier factor (To eventually overlap ghost particles)
    
    Returned value:
    Number of written particles
    """
    n_box = 0
    dist = sep * h
    d = 0.5 * dr
    while d < dist:
        x = r[0] + d * n[0]
        y = r[1] + d * n[1]
        n_box += 1
        write_particle(fhandler, (x, y), n=n, mass=m_fact*refd*dr**2, imove=-1,
                       associated=associated)
        d += dr
    return n_box


def write_wall(fhandler, orig, dest, n, last_id, m_fact=1.0):
    """Write the set of Area boundary elements and Ghost particles along a
    solid wall.
    
    Position arguments:
    orig -- Starting position of the wall line
    dest -- Ending position of the wall line
    n -- Outward normal of the wall
    last_id -- Number of already written particles
    
    Keyword arguments:
    m_fact -- Mass multiplier factor (To eventually overlap ghost particles). It
              is only applied in the extremes of the line (radius sep*h)
    
    Returned value:
    Number of written particles
    """
    n_box = 0
    dist = math.sqrt((dest[0] - orig[0])**2 + (dest[1] - orig[1])**2)
    dir = ((dest[0] - orig[0]) / dist, (dest[1] - orig[1]) / dist)
    d = 0.5 * dr
    while d < dist:
        x = orig[0] + d * dir[0]
        y = orig[1] + d * dir[1]
        imove = -3 if BC == 'BI' else -2
        n_box += 1
        write_particle(fhandler, (x, y), n=n, mass=dr, imove=imove)
        if BC == 'GP':
            associated = last_id + n_box - 1
            mass_factor = 1.0 if sep * h < d < dist - sep*h else m_fact
            n_box += write_gp_particles(fhandler, (x, y), n, associated,
                                        m_fact=mass_factor)
        d += dr
    return n_box

def write_corner(fhandler, r, n, last_id, m_fact=1.0):
    """Write the Ghost particles as well as the common mirror particle, for a
    corner.
    This method is only for convex bodies, and Ghost particles.
    
    Position arguments:
    r -- Mirroring position
    n -- Outward normal of the corner
    last_id -- Number of already written particles
    
    Keyword arguments:
    m_fact -- Mass multiplier factor (To eventually overlap ghost particles)
    
    Returned value:
    Number of written particles
    """
    if BC != 'GP':
        return 0
    n_box = 1
    write_particle(fhandler, r, n=n, mass=dr, imove=-2)
    dir = (math.copysign(1.0, n[0]), math.copysign(1.0, n[1]))
    
    dx = 0.5 * dr
    while dx < sep * h:
        x = r[0] + dir[0] * dx
        dy = 0.5 * dr
        while dy < sep * h:
            y = r[1] + dir[1] * dy
            write_particle(fhandler, (x, y), n=n, mass=m_fact*refd*dr**2,
                           imove=-1, associated=last_id)
            n_box += 1
            dy += dr
        dx += dr
    return n_box

print("Opening fluid particles output file...")
output = open("Fluid.dat", "w")
header = """#############################################################
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
"""
output.write(header)

string = """
    Writing fluid particles...
"""
print(string)

n_fluid = 0
x = 0.5 * dr
while x < L:
    y = 0.5 * dr
    while y < H:
        # Avoid the particles inside the cylinder
        if abs(x - x_square) < 0.5 * D and abs(y - y_square) < 0.5 * D:
            y += dr
            continue
        n_fluid += 1
        write_particle(output, (x, y))
        y += dr
    x += dr
output.close()

string = """
    Writing the outer box...
"""
print(string)
output = open("Box.dat", "w")

n_box = 0
n_box += write_wall(output, (0, 0), (L, 0), ( 0, -1), n_fluid + n_box)
n_box += write_wall(output, (L, 0), (L, H), ( 1,  0), n_fluid + n_box)
n_box += write_wall(output, (L, H), (0, H), ( 0,  1), n_fluid + n_box)
n_box += write_wall(output, (0, H), (0, 0), (-1,  0), n_fluid + n_box)
nmag = 1.0 / math.sqrt(2.0)
n_box += write_corner(output, (0, 0), (-nmag, -nmag), n_fluid + n_box)
n_box += write_corner(output, (L, 0), ( nmag, -nmag), n_fluid + n_box)
n_box += write_corner(output, (L, H), ( nmag,  nmag), n_fluid + n_box)
n_box += write_corner(output, (0, H), (-nmag,  nmag), n_fluid + n_box)
output.close()

string = """
    Writing the inner square...
"""
print(string)
output = open("Square.dat", "w")

n_square = 0
n_square += write_wall(output,
                       (x_square - 0.5 * D, y_square - 0.5 * D),
                       (x_square + 0.5 * D, y_square - 0.5 * D),
                       (0, 1), n_fluid + n_box + n_square, m_fact=0.5)
n_square += write_wall(output,
                       (x_square + 0.5 * D, y_square - 0.5 * D),
                       (x_square + 0.5 * D, y_square + 0.5 * D),
                       (-1, 0), n_fluid + n_box + n_square, m_fact=0.5)
n_square += write_wall(output,
                       (x_square + 0.5 * D, y_square + 0.5 * D),
                       (x_square - 0.5 * D, y_square + 0.5 * D),
                       (0, -1), n_fluid + n_box + n_square, m_fact=0.5)
n_square += write_wall(output,
                       (x_square - 0.5 * D, y_square + 0.5 * D),
                       (x_square - 0.5 * D, y_square - 0.5 * D),
                       (1, 0), n_fluid + n_box + n_square, m_fact=0.5)
output.close()


# It is faster replacing directly N in the files, and reading it with FastASCII,
# than using ASCII reader with the muParser interpreter
fnames = ['Fluid.dat', 'Box.dat', 'Square.dat']
for fname in fnames:
    f = open(fname, 'r')
    txt = f.read()
    f.close()
    txt = txt.replace('N', str(n_fluid + n_box + n_square))
    f = open(fname, 'w')
    f.write(txt)
    f.close()

# XML definition generation
# =========================

templates_path = path.join('@EXAMPLE_DEST_DIR@', 'templates')
XML = ('Main.xml', 'Settings.xml', 'SPH.xml', 'Time.xml', 'Motion.xml',
       'Fluid.xml', 'Box.xml', 'Square.xml', 'BI.xml', 'GP.xml', 'plot_f.py')

domain_min = str(domain_min).replace('(', '').replace(')', '')
domain_max = str(domain_max).replace('(', '').replace(')', '')

data = {'DR':str(dr), 'HFAC':str(hfac), 'CS':str(cs), 'COURANT':str(courant),
        'DOMAIN_MIN':domain_min, 'DOMAIN_MAX':domain_max, 'REFD':str(refd),
        'VISC_DYN':str(visc_dyn), 'DELTA':str(delta), 'G':str(g),
        'NFLUID':str(n_fluid), 'NBOX':str(n_box), 'NSQUARE':str(n_square),
        'NY':str(ny), 'L':str(L), 'H':str(H), 'D':str(D),
        'U':str(U), 'P0':str(p0), 'BC':BC,
        'XSQUARE':str(x_square), 'YSQUARE':str(y_square),}
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