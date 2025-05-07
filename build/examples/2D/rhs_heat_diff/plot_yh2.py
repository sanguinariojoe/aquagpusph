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
import json
import meshio
import numpy as np
#import sodshock
import matplotlib.pyplot as plt
import matplotlib.animation as animation


DR = 0.002
L = 0.2
T = 0.001
GAMMA = 1.4
CS = 349.4923553454445
P = [100000.0, 100000.0]
RHO = [1.156839149618933, 1.1461795999094495]
E = [216106.1026352289, 218115.9043659044]


def read_vtu():
    with open('output.vtu.series') as json_data:
        data = json.load(json_data)
        fname, t = data['files'][-1]['name'], data['files'][-1]['time']
    mesh = meshio.read(fname)
    # We must choose the points.
    mask = np.zeros(len(mesh.points), dtype=bool)
    mask[(np.abs(mesh.points[:, 1]) < 0.25 * DR)] = True
    x = mesh.points[mask, 0]
    rho = mesh.point_data['rho'][mask]
    p = mesh.point_data['p'][mask]
    eint = mesh.point_data['eint'][mask]
    y_H2 = mesh.point_data['y_H2'][mask]
    u = np.linalg.norm(mesh.point_data['u'][mask, :], axis=1)
    sorter = np.argsort(x)
    return t, x[sorter], rho[sorter], p[sorter], eint[sorter], u[sorter], y_H2[sorter]


fig = plt.figure()
ax = fig.add_subplot(111)

# Create the lines
exp, = ax.plot([0.0], [0.0], color="black", linewidth=1.0, linestyle='-')
sph, = ax.plot([0.0], [0.0], color="red", linewidth=1.0, linestyle='--')
# Set some options
ax.grid()
ax.set_xlim(-0.5, 0.5)
ax.set_ylim(0, 1.0)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$x / R$")
ax.set_ylabel(r"$e / e_0$")
ax.set_title(r"$t = 0$")

# Animate
def update(frame_index):
    plt.tight_layout()
    try:
        # SPH data
        t, x, _, _, e, _, y_H2 = read_vtu()
        # Analytic solution
        npts = len(x)
        left_state = (P[0], RHO[0], 0)
        right_state = (P[1], RHO[1], 0.)
        #_, _, exp_data = sodshock.solve(left_state=left_state,
        #                                right_state=right_state,
        #                                geometry=(-0.5 * L, 0.5 * L, 0),
        #                                t=t, 
        #                                gamma=GAMMA,
        #                                npts=npts,
        #                                dustFrac=0.0)
    except IndexError:
        return
    except FileNotFoundError:
        return
    sph.set_data(x / L, y_H2 / max(y_H2))
    #exp.set_data(exp_data['x'] / L, exp_data['energy'] / max(E))
    ax.set_title(r"$t \,\, c_0 / L = {}$".format(t / T))

update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
