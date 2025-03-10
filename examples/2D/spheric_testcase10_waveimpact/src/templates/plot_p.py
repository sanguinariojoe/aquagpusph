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
from scipy.signal import savgol_filter
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


line = None


def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('sensors.out')
        t = data[0]
        p = [d * s * 2 for d, s in zip(data[1], data[2])]
    except IndexError:
        return
    except FileNotFoundError:
        return
    try:
        line.set_data(t, savgol_filter(p, 5, 3))
    except ValueError:
        # Not enough data yet
        line.set_data(t, p)


fig = plt.figure()
ax = fig.add_subplot(111)

FNAME = path.join('@EXAMPLE_DEST_DIR@', 'lateral_water_1x.txt')
T,P,A,DADT,_,_ = np.loadtxt(FNAME,
                            delimiter='\t',
                            skiprows=1,
                            unpack=True)
exp_t = T
exp_p = 100.0 * P

ax.plot(exp_t,
        exp_p,
        label=r'$\mathrm{Experiments}$',
        color="red",
        linewidth=1.0)
t = [0.0]
p = [0.0]
line, = ax.plot(t,
                p,
                label=r'$\mathrm{SPH}$',
                color="black",
                linewidth=1.0)
# Set some options
ax.grid()
ax.legend(loc='best')
ax.set_xlim(0, 5)
ax.set_ylim(-1000, 5000)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t \, [\mathrm{s}]$")
ax.set_ylabel(r"$p \, [\mathrm{Pa}]$")

update(0)
ani = animation.FuncAnimation(fig, update, interval=5000)
plt.show()
