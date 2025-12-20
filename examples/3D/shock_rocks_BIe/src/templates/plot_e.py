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


M = 3900
I = (5560, 5270, 4180)
G = {{G}}

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
ax = fig.add_subplot(111)

line_ek, = ax.plot([0], [0],
                   label=r'$E_k$',
                   color="red",
                   linewidth=1.0,
                   linestyle='-')
line_Ek, = ax.plot([0], [0],
                   label=None,
                   color="red",
                   linewidth=1.0,
                   linestyle='--')
line_ec, = ax.plot([0], [0],
                   label=r'$E_c$',
                   color="blue",
                   linewidth=1.0,
                   linestyle='-')
line_ep, = ax.plot([0], [0],
                   label=r'$E_p$',
                   color="darkgreen",
                   linewidth=1.0,
                   linestyle='-')
line_Ep, = ax.plot([0], [0],
                   label=None,
                   color="darkgreen",
                   linewidth=1.0,
                   linestyle='--')
line_et, = ax.plot([0], [0],
                   label=r'$E_t$',
                   color="black",
                   linewidth=1.0,
                   linestyle='-')
# Set some options
ax.grid()
ax.legend(loc='best')
ax.set_xlim(0.0, {{T}})
ax.set_ylim(-0.1, 0.1)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$t \, [\mathrm{s}]$")
ax.set_ylabel(r"$\mathcal{E}(t) \, [\mathrm{J}]$")


# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        data = readFile('Energy.dat')
        t = data[0]
        ek, ep, ec = data[1], data[2], data[3]
        data = readFile('motion.out')
        z = data[3][:len(t)]
        vx, vy, vz = data[5][:len(t)], data[6][:len(t)], data[7][:len(t)]
        mv2 = [M * (vvx**2 + vvy**2 + vvz**2) \
            for vvx, vvy, vvz in zip(vx, vy, vz)]
        wx, wy, wz = data[17][:len(t)], data[18][:len(t)], data[19][:len(t)]
        iw2 = [(I[0] * wwx**2 + I[1] * wwy**2 + I[2] * wwz**2) \
            for wwx, wwy, wwz in zip(wx, wy, wz)]
        Ep = [M * G * zz for zz in z]
        Ek = [0.5 * (mvv + iww) for mvv, iww in zip(mv2, iw2)]
        Ep = [eep - Ep[0] for eep in Ep]
        Ek = [eek - Ek[0] for eek in Ek]
        et = [ek[i] + ep[i] + ec[i] + Ek[i] + Ep[i] for i in range(len(ek))]
    except IndexError:
        return
    except FileNotFoundError:
        return
    line_ek.set_data(t, ek)
    line_Ek.set_data(t, Ek)
    line_ec.set_data(t, ec)
    line_ep.set_data(t, ep)
    line_Ep.set_data(t, Ep)
    line_et.set_data(t, et)
    ax.set_xlim(0.0, max(t))
    ee = ek + Ek + ec + ep + Ep + et
    emax = max(-min(ee), max(ee))
    ax.set_ylim(min(ee) - 0.05 * emax, max(ee) + 0.05 * emax)


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
