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


data = readFile('britan_etal_1995_x.csv')
t = [t * 1e-3 for t in data[0]]
x = [x * 1e-3 for x in data[1]]
fig, ax = plt.subplots(1, 1)
rx, = ax.plot([0.0], [0.0], color="k", linestyle="-", linewidth=2.0)
ax.plot(t, x, color="r", marker='x', linestyle="--", linewidth=1.0)
ax.grid()
ax.set_xlim(0.0, {{T}})
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t \, [\mathrm{s}]$")
ax.set_ylabel(r"$x \, [\mathrm{m}]$")

# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('motion.out')
        t = data[0]
        x = data[1]
    except IndexError:
        return
    except FileNotFoundError:
        return
    rx.set_data(t, x)
    ax.set_xlim(0, max(t))
    ax.set_ylim(0, max(x))


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
