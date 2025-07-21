# ******************************************************************************
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
# ******************************************************************************
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
# ******************************************************************************

import sys
import os
import json
import meshio
import numpy as np
# import sodshock
import matplotlib.pyplot as plt
import matplotlib.animation as animation


DR = {{DR}}
L = {{L}}
T = {{T}}
GAMMA = {{GAMMA}}
CS = {{CS}}
P = [{{P1}}, {{P2}}]
RHO = [{{RHO1}}, {{RHO2}}]
E = [{{E1}}, {{E2}}]

DT={{DT}}

#DT *= 1012.0

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
    Temp = mesh.point_data['T'][mask]
    u = np.linalg.norm(mesh.point_data['u'][mask, :], axis=1)
    sorter = np.argsort(x)
    return t, x[sorter], rho[sorter], p[sorter], eint[sorter], u[sorter], Temp[sorter],


fig = plt.figure()
ax = fig.add_subplot(111)

# Create the lines
exp, = ax.plot([0.0], [0.0], color="black", linewidth=1.0, linestyle='-')
sph, = ax.plot([0.0], [0.0], color="red", linewidth=1.0, linestyle='--')
# Set some options
ax.grid()
ax.set_xlim(-0.5, 0.5)
#ax.set_ylim(0.9950, 1.0001)
ax.set_autoscale_on(False)
ax.set_xlabel(r"$x / R$")
ax.set_ylabel(r"$T$")
ax.set_title(r"$t = 0$")

# Animate

def solutionferrc(minx, maxx, time, npts):
    
    from scipy import special
    from numpy import linspace, sqrt
    
    x = linspace(minx, maxx, npts)
    
    Temperature = special.erf(x/(2*sqrt(DT * time))) 
    
    return Temperature
    
    
def update(frame_index):
    plt.tight_layout()
    try:
        # SPH data
        t, x, _, _, e, _, T_sph = read_vtu()
        # Analytic solution
        npts = len(x)
        left_state = (P[0], RHO[0], 0)
        right_state = (P[1], RHO[1], 0.)
        # _, _, exp_data = sodshock.solve(left_state=left_state,
        #                                right_state=right_state,
        #                                geometry=(-0.5 * L, 0.5 * L, 0),
        #                                t=t,
        #                                gamma=GAMMA,
        #                                npts=npts,
        #                                dustFrac=0.0)
        
        T_exp = solutionferrc(min(x), max(x), t, npts)
    except IndexError:
        return
    except FileNotFoundError:
        return
    #sph.set_data(x / L, e / max(E))
    # exp.set_data(exp_data['x'] / L, exp_data['energy'] / max(E))
    
    maxT=max(T_sph)
    minT=min(T_sph)
    
    
    
    
    sph.set_data(x / L, T_sph)
    exp.set_data(x / L, minT + (maxT-minT)/2.0 * (T_exp+1.0))
    
    #exp.set_data(x / L, T_exp)
    ax.set_title(r"$t \,\, c_0 / L = {}$".format(t / T))
    ax.set_ylim(minT*0.95, maxT*1.05)


update(0)
ani = animation.FuncAnimation(fig, update, interval=1000)
plt.show()
