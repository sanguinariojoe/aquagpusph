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

import math
import numpy as np
import sys
import os
from os import path
try:
    from PyQt4 import QtGui
except:
    try:
        from PySide import QtGui
    except:
        raise ImportError("PyQt4 or PySide is required to use this tool")

try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
except:
    raise ImportError("matplotlib is required to use this tool")


class FigureController(FigureCanvas):
    """Matplotlib figure widget controller"""

    def __init__(self):
        """Constructor"""
        # Create the figure in the canvas
        self.fig = Figure()
        self.ax_em = self.fig.add_subplot(311)
        self.ax_ec = self.fig.add_subplot(312)
        self.ax_es = self.fig.add_subplot(313)
        FigureCanvas.__init__(self, self.fig)
        # generates first "empty" plot
        self.em, = self.ax_em.plot([0.0],
                                   [0.0],
                                   label=r'$\frac{\mathrm{d} \mathcal{E}_m}{\mathrm{d} t}$',
                                   color="black",
                                   linestyle="-",
                                   linewidth=1.0)
        self.ev, = self.ax_em.plot([0.0],
                                   [0.0],
                                   label=r'$\frac{\mathrm{d} \mathcal{E}_\mu}{\mathrm{d} t}$',
                                   color="red",
                                   linestyle="-",
                                   linewidth=1.0)
        self.ec, = self.ax_ec.plot([0.0],
                                   [0.0],
                                   label=r'$\frac{\mathrm{d} \mathcal{E}_c}{\mathrm{d} t}$',
                                   color="black",
                                   linestyle="-",
                                   linewidth=1.0)
        self.ed, = self.ax_ec.plot([0.0],
                                   [0.0],
                                   label=r'$\frac{\mathrm{d} \mathcal{E}_\delta}{\mathrm{d} t}$',
                                   color="red",
                                   linestyle="-",
                                   linewidth=1.0)
        self.es, = self.ax_es.plot([0.0],
                                   [0.0],
                                   label=r'$\frac{\mathrm{d} \mathcal{E}_{\partial \Omega}}{\mathrm{d} t} - W_{\Omega \rightarrow \partial \Omega}$',
                                   color="red",
                                   linestyle="-",
                                   linewidth=1.0)
        # Set some options
        for ax in (self.ax_em, self.ax_ec, self.ax_es):
            ax.grid()
            ax.legend(loc='lower left')
            ax.set_xlim(0, 1.0)
            ax.set_ylim(-1.0, 1.0)
            ax.set_autoscale_on(False)
            ax.set_xlabel(r"$t$", fontsize=21)
            ax.set_ylabel(r"$\frac{\mathrm{d} \mathcal{E}}{\mathrm{d} t}$", fontsize=21)
        # force the figure redraw
        self.fig.canvas.draw()
        # call the update method (to speed-up visualization)
        self.timerEvent(None)
        # start timer, trigger event every 1000 millisecs (=1sec)
        self.timer = self.startTimer(1000)

    def readFile(self, filepath):
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
                data.append(map(float, fields))
            except:
                continue
        # Transpose the data
        return map(list, zip(*data))

    def timerEvent(self, evt):
        """Custom timerEvent code, called at timer event receive"""
        # Read and plot the new data
        data = self.readFile('Power.dat')
        t = np.asarray(data[0])
        ek = np.asarray(data[2])
        ep = np.asarray(data[3])
        ec = np.asarray(data[4])
        ev = np.asarray(data[5])
        ed = np.asarray(data[6])
        es = np.asarray(data[7])

        data = self.readFile('PressureForces.dat')
        tp = np.asarray(data[0])
        fpx = (np.asarray(data[1]) + np.asarray(data[7]))
        data = self.readFile('ViscousForces.dat')
        tv = np.asarray(data[0])
        fvx = (np.asarray(data[1]) + np.asarray(data[7]))

        data = self.readFile('Motion.dat')
        tm = np.asarray(data[0])
        u = np.asarray(data[2])

        l = min((len(t), len(tp), len(tv), len(tm)))
        t = t[:l]
        em = ek[:l] + ep[:l]
        ec = ec[:l]
        ev = ev[:l]
        ed = ed[:l]
        es = es[:l] - (fpx[:l] + fvx[:l]) * u[:l]

        for ax in (self.ax_em, self.ax_ec, self.ax_es):
            ax.set_xlim(0, max(t))

        em_max = np.max(np.abs(np.concatenate((em, ev), axis=0)))
        self.em.set_data(t, em)
        self.ev.set_data(t, ev)
        self.ax_em.set_ylim(-1.05 * em_max, 1.05 * em_max)

        ec_max = np.max(np.abs(np.concatenate((ec, ed), axis=0)))
        self.ec.set_data(t, ec)
        self.ed.set_data(t, ed)
        self.ax_ec.set_ylim(-1.05 * ec_max, 1.05 * ec_max)

        es_max = np.max(np.abs(es))
        self.es.set_data(t, es)
        self.ax_es.set_ylim(-1.05 * es_max, 1.05 * es_max)

        # Redraw
        self.fig.canvas.draw()


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    widget = FigureController()
    widget.setWindowTitle("Power components")
    widget.show()
    sys.exit(app.exec_())
