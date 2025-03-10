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

t = [0.0]
ek = [0.0]
ec = [0.0]
et = [0.0]
line_ek, = ax.plot(t,
                   ek,
                   label=r'$E_k$',
                   color="red",
                   linewidth=1.0)
line_ep, = ax.plot(t,
                   ec,
                   label=r'$E_p$',
                   color="green",
                   linewidth=1.0)
line_ec, = ax.plot(t,
                   ec,
                   label=r'$E_c$',
                   color="blue",
                   linewidth=1.0)
line_et, = ax.plot(t,
                   et,
                   label=r'$E_t$',
                   color="black",
                   linewidth=1.0)
# Set some options
ax.grid()
ax.legend(loc='best')
ax.set_xlim(0.0, 0.1)
ax.set_ylim(-2.0, 25.0)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$\mathcal{E}(t)$")


# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('Energy.dat')
        t = data[0]
        ek = data[1]
        ep = data[2]
        ec = data[3]
        et = [ek[i] + ep[i] + ec[i] for i in range(len(ek))]
    except IndexError:
        return
    except FileNotFoundError:
        return
    ax.set_xlim(0, t[-1])
    bounds = [min(ek + ep + ec + et), max(ek + ep + ec + et)]
    dbound = 0.05 * max(abs(bounds[0]), abs(bounds[1]))
    ax.set_ylim(bounds[0] - dbound, bounds[1] + dbound)
    line_ek.set_data(t, ek)
    line_ep.set_data(t, ep)
    line_ec.set_data(t, ec)
    line_et.set_data(t, et)


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
