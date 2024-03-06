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


line = None


fig = plt.figure()
ax = fig.add_subplot(111)

# Load the percentiles
data = readFile('Fig21_h30_sensor1_pressure_th.csv')
ax.plot(data[0], data[1], label=None, color="black", linewidth=1.0,
        linestyle='--')
ax.plot(data[0], data[3], label=None, color="black", linewidth=1.0,
        linestyle='--')
# Create the SPH line
t = [0.0]
p = [0.0]
line, = ax.plot(t, p, label=r'$p(t)$', color="red", linewidth=1.0)
# Set some options
ax.grid()
T = ({{H}} / {{G}})**0.5
P = {{REFD}} * {{G}} * {{H}}
ax.set_xlim(0.0, {{T}} / T)
ax.set_ylim(-0.3, 2.5)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t / T$")
ax.set_ylabel(r"$p(t) / \rho_0 g H$")

# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('sensors.out')
        t = [d / T for d in data[0]]
        p = [d / P for d in data[2]]
    except IndexError:
        return
    except FileNotFoundError:
        return
    line.set_data(t, p)

update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
