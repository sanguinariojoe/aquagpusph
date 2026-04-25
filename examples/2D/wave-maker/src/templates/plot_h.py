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

import os
from os import path
import numpy as np
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
        l = l.strip()
        while l.find('  ') != -1:
            l = l.replace('  ', ' ')
        fields = l.split(' ')
        try:
            data.append(map(float, fields))
        except:
            continue
    # Transpose the data
    return [list(d) for d in zip(*data)]


lines = []


titles = ('H1', 'H2', 'H3', 'H4')
fig = plt.figure()
ax11 = fig.add_subplot(221)
ax21 = fig.add_subplot(222, sharey=ax11)
ax12 = fig.add_subplot(223, sharex=ax11)
ax22 = fig.add_subplot(224, sharex=ax21, sharey=ax12)
axes = (ax11, ax21, ax12, ax22)


def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('sensors_h.out')
        t = data[0]
        hh = (data[-4], data[-3], data[-2], data[-1])
    except IndexError:
        return
    except FileNotFoundError:
        return
    min_h, max_h = np.min(hh), np.max(hh)
    dh = max_h - min_h
    for i, h in enumerate(hh):
        lines[i].set_data(t, h)
        axes[i].set_xlim(0, max(t))
        axes[i].set_ylim(min_h - 0.05 * dh, max_h + 0.05 * dh)


for i, ax in enumerate(axes):
    t = [0.0]
    h = [0.0]
    line, = ax.plot(t, h, color="black", linewidth=1.0)
    lines.append(line)
    # Set some options
    ax.grid()
    ax.set_title(titles[i])
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, 0.1)
    ax.set_autoscale_on(False)
    if i > 1:
        ax.set_xlabel(r"$t \, [\mathrm{s}]$")
    else:
        plt.setp(ax.get_xticklabels(), visible=False)
    if i in (0, 2):
        ax.set_ylabel(r"$H \, [\mathrm{m}]$")
    else:
        plt.setp(ax.get_yticklabels(), visible=False)

update(0)
ani = animation.FuncAnimation(fig, update, interval=5000)
plt.show()
