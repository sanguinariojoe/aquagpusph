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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# Some case constants
rho = {{REFD}}
D = {{D}}
U = {{U}}
# Coefficients adimensionalization factor
COEFF_FAC = 1.0 / (0.5 * rho * D * U**2)
TIME_FAC = 1


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


fig = plt.figure()
data = readFile('Force_Re100.dat')
t = np.asarray(data[0])
fx = np.asarray(data[1]) + np.asarray(data[2])
exp, = plt.plot(t, fx,
                 label=r'$\mathrm{FDM}$',
                 color="red",
                 linestyle="-",
                 linewidth=1.0)
cd, = plt.plot([0], [0],
               label=r'$\mathrm{AQUAgpusph}$',
               color="black",
               linestyle="-",
               linewidth=1.0)
ax = plt.gca()
ax.grid()
ax.legend(loc='upper right')
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$c_D$")


# Animate
def update(frame_index):
    try:
        data = readFile('PressureForces.dat')
        tp = np.asarray(data[0]) * TIME_FAC
        fpx = (np.asarray(data[1]) + np.asarray(data[7])) * COEFF_FAC
        data = readFile('ViscousForces.dat')
        tv = np.asarray(data[0]) * TIME_FAC
        fvx = (np.asarray(data[1]) + np.asarray(data[7])) * COEFF_FAC
        t = tp if len(tp) < len(tv) else tv
        fx = fpx[:len(t)] + fvx[:len(t)]

        cd.set_data(t, -fx)
        cdmax = np.max(np.abs(fx))
        ax = plt.gca()
        ax.set_xlim(0, np.max(t))
        ax.set_ylim(-1.05 * cdmax, 1.05 * cdmax)
        plt.tight_layout()
    except IndexError:
        return
    except FileNotFoundError:
        return


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
