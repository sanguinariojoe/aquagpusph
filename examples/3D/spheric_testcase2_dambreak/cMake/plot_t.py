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

import math
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


fig = plt.figure()
ax = fig.add_subplot(111)

t = [0.0]
e = [0.0]
fove = ax.fill_between(t,
                       0,
                       e,
                       facecolor='red',
                       linewidth=0.0)
fave = ax.fill_between(t,
                       0,
                       e,
                       facecolor='blue',
                       linestyle="-",
                       linewidth=0.0)
love, = ax.plot(t,
                       e,
                       color='#990000',
                       linestyle="-",
                       linewidth=2.0,
                       label='Average overhead')
lave, = ax.plot(t,
                       e,
                       color='#000099',
                       linestyle="-",
                       linewidth=2.0,
                       label='Average elapsed')
line, = ax.plot(t,
                       e,
                       color="black",
                       linestyle="-",
                       linewidth=1.0,
                       alpha=0.5,
                       label='Elapsed')


def update(frame_index):
    plt.tight_layout()
    data = readFile('Performance.dat')
    t = data[0]
    e = data[1]
    e_ela = data[2]
    e_ove = data[5]
    # Clear nan values
    for i in range(len(e_ela)):
        if math.isnan(e_ela[i]):
            e_ela[i] = 0.0
        if math.isnan(e_ove[i]):
            e_ove[i] = 0.0
    e_ave = [e_ela[i] - e_ove[i] for i in range(len(e_ela))]
    # clear the fills
    for coll in (ax.collections):
        ax.collections.remove(coll)
    fove = ax.fill_between(t,
                           0,
                           e_ela,
                           facecolor='red',
                           linestyle="-",
                           linewidth=2.0)
    fave = ax.fill_between(t,
                           0,
                           e_ave,
                           facecolor='blue',
                           linestyle="-",
                           linewidth=2.0)
    love.set_data(t, e_ela)
    lave.set_data(t, e_ave)
    line.set_data(t, e)

    ax.set_xlim(0, t[-1])
    ax.set_ylim(0, 1.5 * e_ela[-1])


# Set some options
ax.grid()
ax.set_xlim(0, 0.1)
ax.set_ylim(-0.1, 0.1)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t \, [\mathrm{s}]$", fontsize=21)
ax.set_ylabel(r"$t_{CPU} \, [\mathrm{s}]$", fontsize=21)
ax.legend(handles=[lave, love, line], loc='upper right')

ani = animation.FuncAnimation(fig, update, interval=5000)
plt.show()
