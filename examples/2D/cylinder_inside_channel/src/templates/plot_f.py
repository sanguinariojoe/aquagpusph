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


rho = {{REFD}}
D = {{D}}
U = {{U}}
Re = {{RE}}
# Coefficients adimensionalization factor
COEFF_FAC = 1.0 / (0.5 * rho * D * U**2)
TIME_FAC = U / D


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
    for l in lines[:-1]:  # Skip the last line, which may be unready
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
            data.append([float(f) for f in fields])
        except:
            continue
    # Transpose the data
    return [list(d) for d in zip(*data)]


# Find the Cd from the F. M. White experiments
fexp = open('F.M.White.csv', "r")
lines = fexp.readlines()
fexp.close()
re, cd = [], []
for line in lines:
    line = line.strip()
    if line == "":
        continue
    line = line.replace(" ", "")
    fields = [float(field) for field in line.split(",")]
    if fields[0] >= Re:
        re.append(fields[0])
        cd.append(fields[1])
        break
    re, cd = [fields[0]], [fields[1]]
factor = (Re - re[0]) / (re[1] - re[0])
cd_exp = (1 - factor) * cd[0] + factor * cd[1]

fig = plt.figure()
ax_cl = fig.add_subplot(111)
ax_cd = ax_cl.twinx()
line_cl, = ax_cl.plot([0], [0], color="black", linewidth=1.0)
line_cd, = ax_cd.plot([0], [0], color="red", linewidth=1.0)
line_cd_exp, = ax_cd.plot([0], [0], color="red", linestyle='--', linewidth=1.0)

# Set some options
ax_cl.grid()
ax_cl.set_xlim(0, 5)
ax_cl.set_ylim(-1000, 5000)
ax_cl.set_autoscale_on(False)
ax_cl.set_xlabel(r"$t \, [\mathrm{s}]$")
ax_cl.set_ylabel(r"$C_L$")
ax_cd.set_xlim(0, 5)
ax_cd.set_ylim(-1000, 5000)
ax_cd.set_autoscale_on(False)
ax_cd.set_ylabel(r"$C_D$")

def update(frame_index):
    plt.tight_layout()
    data = readFile('force_iset.out')
    t = data[0]
    fpx = data[1]
    fpy = data[2]
    fvx = data[7]
    fvy = data[8]
    fx = []
    fy = []
    for i in range(len(t)):
        t[i] *= TIME_FAC
        fx.append((fpx[i] + fvx[i]) * COEFF_FAC)
        fy.append((fpy[i] + fvy[i]) * COEFF_FAC)
    line_cl.set_data(t, fy)
    line_cd.set_data(t, fx)
    line_cd_exp.set_data([t[0], t[-1]], [cd_exp, cd_exp])

    clmax = max(max(fy), abs(min(fy)))
    cdmax = max(max(fx), abs(min(fx)))
    ax_cl.set_xlim(0, 1.05 * max(t))
    ax_cd.set_xlim(0, 1.05 * max(t))
    ax_cl.set_ylim(-1.05 * clmax, 1.05 * clmax)
    ax_cd.set_ylim(-1.05 * cdmax, 1.05 * cdmax)
    for tl in ax_cd.get_yticklabels():
        tl.set_color("red")
        # Redraw
        fig.canvas.draw()

update(0)
ani = animation.FuncAnimation(fig, update, interval=5000)
plt.show()
