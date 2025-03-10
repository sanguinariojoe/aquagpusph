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
import quasi_static


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
fp = [0.0]
fe = [0.0]
ft = [0.0]
line_fp, = ax.plot(t,
                   fp,
                   label=r'$\boldsymbol{f}_p$',
                   color="red",
                   linewidth=1.0)
line_fe, = ax.plot(t,
                   fe,
                   label=r'$\boldsymbol{f}_{elastic}$',
                   color="blue",
                   linewidth=1.0)
line_ft, = ax.plot(t,
                   ft,
                   label=r'$\boldsymbol{f}_t$',
                   color="black",
                   linewidth=1.0)
t, _, f, _ = quasi_static.simulate()
t = [d / {{T0}} for d in t]
f = [d / {{F}} for d in f]
ax.plot(t, f, color="black", linewidth=1.0, linestyle='--')
# Set some options
ax.grid()
ax.legend(loc='best')
ax.set_xlim(0.0, {{T}} / {{T0}})
ax.set_ylim(0.0, 2.1)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t / T$")
ax.set_ylabel(r"$f(t) / F$")


# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('force.out')
        t = [d / {{T0}} for d in data[0]]
        fp = [d / {{F}} for d in data[1]]
        fe = [d / {{F}} for d in data[3]]
        ft = [fp[i] + fe[i] for i in range(len(fp))]
    except IndexError:
        return
    except FileNotFoundError:
        return
    line_fp.set_data(t, fp)
    line_fe.set_data(t, fe)
    line_ft.set_data(t, ft)


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
