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
ax_em = fig.add_subplot(311)
ax_ec = fig.add_subplot(312)
ax_es = fig.add_subplot(313)
em, = ax_em.plot([0], [0],
                 label=r'$\frac{\mathrm{d} \mathcal{E}_m}{\mathrm{d} t}$',
                 color="black",
                 linestyle="-",
                 linewidth=1.0)
ev, = ax_em.plot([0], [0],
                 label=r'$\frac{\mathrm{d} \mathcal{E}_\mu}{\mathrm{d} t}$',
                 color="red",
                 linestyle="-",
                 linewidth=1.0)
ec, = ax_ec.plot([0], [0],
                 label=r'$\frac{\mathrm{d} \mathcal{E}_c}{\mathrm{d} t}$',
                 color="black",
                 linestyle="-",
                 linewidth=1.0)
ed, = ax_ec.plot([0], [0],
                 label=r'$\frac{\mathrm{d} \mathcal{E}_\delta}{\mathrm{d} t}$',
                 color="red",
                 linestyle="-",
                 linewidth=1.0)
es, = ax_es.plot([0], [0],
                 label=r'$\frac{\mathrm{d} \mathcal{E}_{\partial \Omega}}{\mathrm{d} t} - W_{\Omega \rightarrow \partial \Omega}$',
                 color="red",
                 linestyle="-",
                 linewidth=1.0)
# Set some options
for ax in (ax_em, ax_ec, ax_es):
    ax.grid()
    ax.legend(loc='lower left')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_autoscale_on(False)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\frac{\mathrm{d} \mathcal{E}}{\mathrm{d} t}$")


# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('Power.dat')
        t = np.asarray(data[0])
        data_ek = np.asarray(data[2])
        data_ep = np.asarray(data[3])
        data_ec = np.asarray(data[4])
        data_ev = np.asarray(data[5])
        data_ed = np.asarray(data[6])
        data_es = np.asarray(data[7])

        data = readFile('PressureForces.dat')
        tp = np.asarray(data[0])
        fpx = (np.asarray(data[1]) + np.asarray(data[7]))
        data = readFile('ViscousForces.dat')
        tv = np.asarray(data[0])
        fvx = (np.asarray(data[1]) + np.asarray(data[7]))

        data = readFile('Motion.dat')
        tm = np.asarray(data[0])
        u = np.asarray(data[2])
    except IndexError:
        return
    except FileNotFoundError:
        return

    l = min((len(t), len(tp), len(tv), len(tm)))
    t = t[:l]
    data_em = data_ek[:l] + data_ep[:l]
    data_ec = data_ec[:l]
    data_ev = data_ev[:l]
    data_ed = data_ed[:l]
    data_es = data_es[:l] - (fpx[:l] + fvx[:l]) * u[:l]

    for ax in (ax_em, ax_ec, ax_es):
        ax.set_xlim(0, max(t))

    em_max = np.max(np.abs(np.concatenate((data_em, data_ev), axis=0)))
    em.set_data(t, data_em)
    ev.set_data(t, data_ev)
    ax_em.set_ylim(-1.05 * em_max, 1.05 * em_max)

    ec_max = np.max(np.abs(np.concatenate((data_ec, data_ed), axis=0)))
    ec.set_data(t, data_ec)
    ed.set_data(t, data_ed)
    ax_ec.set_ylim(-1.05 * ec_max, 1.05 * ec_max)

    es_max = np.max(np.abs(data_es))
    es.set_data(t, data_es)
    ax_es.set_ylim(-1.05 * es_max, 1.05 * es_max)


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
