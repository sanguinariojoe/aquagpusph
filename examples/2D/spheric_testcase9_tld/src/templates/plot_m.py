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
        fields = l.split('\t')
        try:
            data.append(map(float, fields))
        except:
            continue
    # Transpose the data
    return list(map(list, zip(*data)))


line = None


fig = plt.figure()
ax = fig.add_subplot(111)

exp_line, = ax.plot([0.0], [0.0],
                    label=r'$\theta_{Exp}$',
                    color="black",
                    linestyle="--",
                    linewidth=1.0)

line, = ax.plot([0.0], [0.0],
                label=r'$\theta_{SPH}$',
                color="black",
                linestyle="-",
                linewidth=1.0)
# Set some options
ax.grid()
ax.legend(loc='upper left')
ax.set_xlim(0.0, 0.1)
ax.set_ylim(-0.1, 0.1)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t \, [\mathrm{s}]$")
ax.set_ylabel(r"$\theta \, [\mathrm{deg}]$")

# Create a secondary set of axes for the moment
ax2 = ax.twinx()
mline, = ax2.plot([0.0], [0.0],
                  label=r'$M_{SPH}$',
                  color="blue",
                  linewidth=1.0)
# Set some options
ax2.set_xlim(0, 0.1)
ax2.set_ylim(-0.1, 0.1)
ax2.set_autoscale_on(False)
ax2.set_ylabel(r"$M_{fluid} \, [\mathrm{N} \cdot \mathrm{m}]$", color="blue")
for tl in ax2.get_yticklabels():
    tl.set_color("blue")


# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('Motion.dat')
        t = data[0]
        exp_a = data[2]
        a = data[3]
        m = data[5]
    except IndexError:
        return
    except FileNotFoundError:
        return
    exp_line.set_data(t, exp_a)
    line.set_data(t, a)
    mline.set_data(t, m)

    ax.set_xlim(0, t[-1])
    ymax_exp = max(abs(max(exp_a)), abs(min(exp_a)))
    ymax_sph = max(abs(max(a)), abs(min(a)))
    ymax = max((ymax_exp, ymax_sph))
    ax.set_ylim(-1.1 * ymax, 1.1 * ymax)

    ax2.set_xlim(0, t[-1])
    ymax = max(abs(max(m)), abs(min(m)))
    ax2.set_ylim(-1.1 * ymax, 1.1 * ymax)
    for tl in ax2.get_yticklabels():
        tl.set_color("blue")


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
