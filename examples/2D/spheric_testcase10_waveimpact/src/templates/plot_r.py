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
import math
from os import path
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


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
    for l in lines[1:-1]:  # Skip the last line, which may be unready
        l = l.replace('\t', ' ')
        l = l.replace('(', ' ')
        l = l.replace(')', ' ')
        l = l.replace(',', ' ')
        while l.find('  ') != -1:
            l = l.replace('  ', ' ')
        fields = l.strip().split(' ')
        if len(data) and len(fields) != len(data[-1]):
            # Probably the system is writing it
            continue
        try:
            data.append(list(map(float, fields)))
        except:
            continue
    # Transpose the data
    return list(map(list, zip(*data)))


def stats(residues):
    """ Compute the minimum, maximum and average value for each iteration
    :param residues The stack of residues grouped by iteration
    """
    y_min = []
    y = []
    y_max = []
    for it in range(len(residues)):
        y_min.append(np.percentile(residues[it], 0.05))
        y_max.append(np.percentile(residues[it], 0.95))
        y.append(np.percentile(residues[it], 0.5))
    return y_min, y, y_max

fig = plt.figure()
ax = fig.add_subplot(111)

# Create the SPH line
t = [0.0]
r = [0.0]
line_min, = ax.plot([0.0], [0.0], color="black", linewidth=1.0, linestyle='--')
line, = ax.plot([0.0], [0.0], color="black", linewidth=1.0, linestyle='-')
line_max, = ax.plot([0.0], [0.0], color="black", linewidth=1.0, linestyle='--')
# Set some options
ax.grid()
ax.set_xlim(0.0, 1.0)
ax.set_ylim(1e-8, 1.0)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$iter$")
ax.set_ylabel(r"$R_{\Delta t}$")
ax.set_yscale('log')

# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('midpoint.out')
        # Stack the values according to its own iteration
        r = []
        for i in range(len(data[0])):
            it = int(data[1][i])
            if len(r) < it + 1:
                r.append([])
            r[it].append(data[-1][i])
        # Get the stats for each iteration
        x = list(range(len(r)))
        y_min, y, y_max = stats(r)
        ax.set_xlim(0.0, x[-1])
        ax.set_ylim(min(y_min), max(y_max))
    except IndexError:
        return
    except FileNotFoundError:
        return
    line_min.set_data(x, y_min)
    line.set_data(x, y)
    line_max.set_data(x, y_max)

update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
