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


lines = []


def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('sensors_0.out')
        t = data[0]
        pp = data[1:]
    except IndexError:
        return
    except FileNotFoundError:
        return
    for i, p in enumerate(pp):
        try:
            lines[i].set_data(t, savgol_filter(p, 71, 3))
        except ValueError:
            # Not enough data yet
            lines[i].set_data(t, p)


fig = plt.figure()
ax11 = fig.add_subplot(241)
ax12 = fig.add_subplot(242, sharey=ax11)
ax13 = fig.add_subplot(243, sharey=ax11)
ax14 = fig.add_subplot(244, sharey=ax11)
ax21 = fig.add_subplot(245, sharex=ax11)
ax22 = fig.add_subplot(246, sharex=ax12, sharey=ax21)
ax23 = fig.add_subplot(247, sharex=ax13, sharey=ax21)
ax24 = fig.add_subplot(248, sharex=ax14, sharey=ax21)
axes = (ax11, ax12, ax13, ax14,
        ax21, ax22, ax23, ax24)

FNAME = path.join('@EXAMPLE_DEST_DIR@', 'test_case_2_exp_data.dat')
T,P1,P2,P3,P4,P5,P6,P7,P8,_,_,_,_, = readFile(FNAME)
exp_t = T
exp_p = (P1, P2, P3, P4, P5, P6, P7, P8)
titles = ('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8')

for i, ax in enumerate(axes):
    ax.plot(exp_t,
            exp_p[i],
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
    lines.append(line)
    # Set some options
    ax.grid()
    ax.legend(loc='best')
    ax.set_title(titles[i])
    ax.set_xlim(0, 6)
    ax.set_autoscale_on(False)
    if i > 3:
        ax.set_ylim(-2000, 4000)
        ax.set_xlabel(r"$t \, [\mathrm{s}]$")
    else:
        ax.set_ylim(-1000, 15000)
        plt.setp(ax.get_xticklabels(), visible=False)
    if i in (0, 4):
        ax.set_ylabel(r"$p \, [\mathrm{Pa}]$")
    else:
        plt.setp(ax.get_yticklabels(), visible=False)

update(0)
ani = animation.FuncAnimation(fig, update, interval=5000)
plt.show()
