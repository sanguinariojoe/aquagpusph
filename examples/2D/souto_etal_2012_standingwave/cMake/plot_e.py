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
    return map(list, zip(*data))


line = None


fig = plt.figure()
ax = fig.add_subplot(111)

data = readFile('Ekin.dat')
T = data[0][-1] / data[1][-1]
ax.plot(data[1],
        data[3],
        label=r'$\mathrm{Theoric}$',
        color="red",
        linewidth=1.0)

t = [0.0]
e = [0.0]
line, = ax.plot(t,
                e,
                label=r'$\mathrm{SPH}$',
                color="black",
                linewidth=1.0)
# Set some options
ax.grid()
ax.legend(loc='best')
ax.set_xlim(0.0, 8.0)
ax.set_ylim(0.0, 1.0)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t / T$")
ax.set_ylabel(r"$\mathcal{E}_{k}(t) / \mathcal{E}_{k}(0)$")


# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('Energy.dat')
        t = data[0]
        e = data[1]
        e0 = e[0]
        for i in range(len(t)):
            t[i] /= T
            e[i] /= e0
    except IndexError:
        return
    except FileNotFoundError:
        return
    line.set_data(t, e)


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
