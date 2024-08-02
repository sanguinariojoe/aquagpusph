#******************************************************************************
#                                                                             *
#              *    **   *  *   *                           *                 *
#             * *  *  *  *  *  * *                          *                 *
#            ***** *  *  *  * *****  **  ***  *  *  ** ***  ***               *
#            *   * *  *  *  * *   * *  * *  * *  * *   *  * *  *              *
#            *   * *  *  *  * *   * *  * *  * *  *   * *  * *  *              *
#            *   *  ** *  **  *   *  *** ***   *** **  ***  *  *              *
#                                      * *             *                      *
#                                    **  *             *                      *
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

import sys
import os
from os import path
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation


G = {{G}}
ACCX, ACCY, ACCZ, ACCA = 5.9, 0.0, 33.5, 168


def readFile(filepath):
    """ Read and extract data from a file
    :param filepath File ot read
    """
    abspath = filepath
    if not path.isabs(filepath):
        abspath = path.join(path.dirname(path.abspath(__file__)), filepath)
    # Read the file by lines
    f = open(abspath, "r")
    lines = f.readlines()
    f.close()
    data = []
    for l in lines[:-1]:  # Skip the last line, which may be unready
        l = l.strip()
        if l.startswith('#') or l == '':
            continue
        l = l.replace('(', ' ').replace(')', ' ').replace(',', ' ')
        while l.find('  ') != -1:
            l = l.replace('  ', ' ')
        l = l.strip()
        fields = l.split(' ')
        try:
            data.append(map(float, fields))
        except:
            continue
    # Transpose the data
    return list(map(list, zip(*data)))


def rotate_pitch(ax, az, theta):
    rx = [x * math.cos(t) - z * math.sin(t) for x,z,t in zip(ax, az, theta)]
    rz = [z * math.cos(t) + x * math.sin(t) for x,z,t in zip(ax, az, theta)]
    return rx, rz


fig, ((ax_accx, ax_accy), (ax_accz, ax_acca)) = plt.subplots(2, 2)
ax_accx.plot([0.0, {{T}}], [ACCX, ACCX],
             color="black",
             linestyle="--",
             linewidth=1.0)
accx, = ax_accx.plot([0.0], [0.0],
                     color="black",
                     linestyle="-",
                     linewidth=1.0)
t, a = readFile('rsc/accy.csv')
ax_accy.plot([0.0, {{T}}], [ACCY, ACCY],
             color="black",
             linestyle="--",
             linewidth=1.0)
accy, = ax_accy.plot([0.0], [0.0],
                     color="black",
                     linestyle="-",
                     linewidth=1.0)
t, a = readFile('rsc/accz.csv')
ax_accz.plot([0.0, {{T}}], [ACCZ, ACCZ],
             color="black",
             linestyle="--",
             linewidth=1.0)
accz, = ax_accz.plot([0.0], [0.0],
                     color="black",
                     linestyle="-",
                     linewidth=1.0)
t, a = readFile('rsc/acca.csv')
ax_acca.plot([0.0, {{T}}], [ACCA, ACCA],
             color="black",
             linestyle="--",
             linewidth=1.0)
acca, = ax_acca.plot([0.0], [0.0],
                     color="black",
                     linestyle="-",
                     linewidth=1.0)
# Set some options
for ax in (ax_accx, ax_accy, ax_accz, ax_acca):
    ax.grid()
    ax.set_xlim(0.0, {{T}})
    ax.set_autoscale_on(False)
    ax.set_xlabel(r"$t \, [\mathrm{s}]$")
ax_accx.set_ylabel(r"$\ddot{x} \,\, \vert \boldsymbol{g} \vert^{-1}$")
ax_accy.set_ylabel(r"$\ddot{y} \,\, \vert \boldsymbol{g} \vert^{-1}$")
ax_accz.set_ylabel(r"$\ddot{z} \,\, \vert \boldsymbol{g} \vert^{-1}$")
ax_acca.set_ylabel(r"$\ddot{\theta} \,\, [\mathrm{rad} \cdot \mathrm{s}^{-2}]$")

# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('motion.out')
        t = data[0]
        ax = [a / G for a in data[9]]
        ay = [a / G for a in data[10]]
        az = [a / G for a in data[11]]
        theta = data[14]
        aa = data[22]
        ax, az = rotate_pitch(ax, az, theta)
    except IndexError:
        return
    except FileNotFoundError:
        return
    accx.set_data(t, ax)
    ax_accx.set_xlim(0, max(t))
    ax_accx.set_ylim(min(ax) - 0.05 * ACCX, max(max(ax), ACCX) + 0.05 * ACCX)
    accy.set_data(t, ay)
    ax_accy.set_xlim(0, max(t))
    ax_accy.set_ylim(min(ay) - 0.05 * ACCY, max(max(ay), ACCY) + 0.05 * ACCY)
    accz.set_data(t, az)
    ax_accz.set_xlim(0, max(t))
    ax_accz.set_ylim(min(az) - 0.05 * ACCZ, max(max(az), ACCZ) + 0.05 * ACCZ)
    acca.set_data(t, aa)
    ax_acca.set_xlim(0, max(t))
    ax_acca.set_ylim(min(aa) - 0.05 * ACCA, max(max(aa), ACCA) + 0.05 * ACCA)


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
