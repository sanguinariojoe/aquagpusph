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
ek = [0.0]
ec = [0.0]
ew = [0.0]
et = [0.0]
line_ek, = ax.plot(t,
                   ek,
                   label=r'$E_k$',
                   color="red",
                   linewidth=1.0)
line_ec, = ax.plot(t,
                   ec,
                   label=r'$E_c$',
                   color="blue",
                   linewidth=1.0)
line_ew, = ax.plot(t,
                   ew,
                   label=r'$W$',
                   color="green",
                   linewidth=1.0)
line_et, = ax.plot(t,
                   et,
                   label=r'$E_t$',
                   color="black",
                   linewidth=1.0)
t, _, _, e = quasi_static.simulate()
t = [d / {{T0}} for d in t]
e = [d / abs({{F}} * {{X0}}) for d in e]
ax.plot(t, e, color="black", linewidth=1.0, linestyle='--')
# Set some options
ax.grid()
ax.legend(loc='best')
ax.set_xlim(0.0, {{T}} / {{T0}})
ax.set_ylim(-2.1, 2.1)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t / T$")
ax.set_ylabel(r"$\mathcal{E}(t) / (F \cdot x(t=0))$")


# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('Energy.dat')
        t = [d / {{T0}} for d in data[0]]
        ek = [d / abs({{F}} * {{X0}}) for d in data[1]]
        ec = [d / abs({{F}} * {{X0}}) for d in data[3]]
        data = readFile('spring.out')
        x = [d - ({{X0}}) for d in data[1]]
        dxdt = data[2]
        ew = [(0.5 * {{M}} * dxdt[i]**2 + {{F}} * x[i]) / abs({{F}} * {{X0}}) \
            for i in range(len(t))]
        et = [ek[i] + ec[i] + ew[i] for i in range(len(ek))]
    except IndexError:
        return
    except FileNotFoundError:
        return
    line_ek.set_data(t, ek)
    line_ec.set_data(t, ec)
    line_ew.set_data(t, ew)
    line_et.set_data(t, et)


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
