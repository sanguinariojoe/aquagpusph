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
from math import *
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


g = 9.81
h = 0.3
l = 0.6
L = 1.61
rho = 997.0


class FigureController(FigureCanvas):
    """Matplotlib figure widget controller"""

    def __init__(self):
        """Constructor"""
        self.ax = [None, None, None, None]
        self.line = [None, None, None, None]
        self.t = [[0.0], [0.0], [0.0], [0.0]]
        self.p = [[0.0], [0.0], [0.0], [0.0]]
        # Create the figure in the canvas
        self.fig = Figure()
        self.ax[0] = self.fig.add_subplot(221)
        self.ax[1] = self.fig.add_subplot(222)
        self.ax[2] = self.fig.add_subplot(223)
        self.ax[3] = self.fig.add_subplot(224)
        FigureCanvas.__init__(self, self.fig)
        # generates first "empty" plot
        files = (('Fig20_filtered_2.dat', 'Fig20_filtered_3.dat'),
                 ('Fig30_filtered_1.dat', 'Fig30_filtered_2.dat'),
                 ('Fig31_filtered_3.dat', 'Fig31_filtered_4.dat'),
                 ('Fig29_filtered_3.dat', 'Fig29_filtered_4.dat'))
        for i,f in enumerate(files):
            data = self.readFile(f[0])
            if i in [0]:  # List of figures that nned to be renormalized
                t = [data[0][j] * (g / h) for j in range(len(data[1]))]
                p = [data[1][j] * 100 for j in range(len(data[1]))]
            else:
                t = data[0]
                p = data[1]
            self.ax[i].plot(t,
                            p,
                            label=r'2.5$\%$ and 97.5$\%$ percentiles',
                            color="black",
                            linestyle='--',
                            linewidth=2.0)
            data = self.readFile(f[1])
            if i in [0]:  # List of figures that nned to be renormalized
                t = [data[0][j] * (g / h) for j in range(len(data[1]))]
                p = [data[1][j] * 100 for j in range(len(data[1]))]
            else:
                t = data[0]
                p = data[1]
            self.ax[i].plot(t,
                            p,
                            label='_nolegend_',
                            color="black",
                            linestyle='--',
                            linewidth=2.0)
            
            self.line[i], = self.ax[i].plot(self.t[i],
                                            self.p[i],
                                            label='SPH',
                                            linewidth=1.0)
            self.ax[i].grid()
            self.ax[i].legend(loc='best')
            self.ax[i].set_autoscale_on(False)
            self.ax[i].set_xlabel(r"$t \sqrt{\frac{g}{H}}$", fontsize=21)
            self.ax[i].set_ylabel(r"$\frac{p}{\rho g H}$", fontsize=21)
            self.ax[i].set_xlim(0, 7)
            self.ax[i].set_ylim(-0.2, 4)
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
        for l in lines:
            l = l.strip()
            fields = l.split('\t')
            try:
                data.append(map(float, fields))
            except:
                continue
        # Transpose the data
        return map(list, zip(*data))

    def timerEvent(self, evt):
        """Custom timerEvent code, called at timer event receive"""
        # Read and plot the new data
        data = self.readFile('../Sensors.dat')
        for i in range(0, 4):
            self.t[i] = [data[0][j] * sqrt(g / h) for j in range(len(data[0]))]
            self.p[i] = []
            for j in range(len(data[0])):
                self.p[i].append(data[3 + 6 * i][j] / (rho * g * h))
                sumW = data[6 + 6 * i][j]
                # Apply a threshold to avoid the interpolation based just on
                # free surface particles
                if sumW < 0.3:
                    self.p[i][j] = 0.0
            self.line[i].set_data(self.t[i], self.p[i])
        # Redraw
        self.fig.canvas.draw()


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    widget = FigureController()
    widget.setWindowTitle("Sensor 4")
    widget.show()
    sys.exit(app.exec_())
