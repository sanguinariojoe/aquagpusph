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
    import matplotlib.pyplot as plt
except:
    raise ImportError("matplotlib is required to use this tool")


cs = 45.0
gamma = 1.0
refd = 998.0


def eos(dens):
    return cs * cs * refd / gamma * ((dens / refd)**gamma - 1.0)


class FigureController(FigureCanvas):
    """Matplotlib figure widget controller"""

    def __init__(self):
        """Constructor"""
        # Create the figure in the canvas
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
        # generates first "empty" plot
        data = self.readFile('lateral_water_1x.txt')
        self.exp_t = data[0]
        self.exp_p = [100.0 * data[1][i] for i in range(len(data[1]))]
        self.exp_line, = self.ax.plot(self.exp_t,
                                      self.exp_p,
                                      label=r'$p_{exp}$',
                                      color="red",
                                      linewidth=1.0)
        self.t = [0.0]
        self.p = [0.0]
        self.pmin = [0.0]
        self.pmax = [0.0]
        self.var_line = self.ax.fill_between(self.t,
                                             self.pmin,
                                             self.pmax,
                                             interpolate=True,
                                             facecolor=(0.75, 0.75, 0.75, 1.0),
                                             edgecolor=(0.5, 0.5, 0.5, 1.0),
                                             label=r'$p_{SPH} \pm \sigma\left(p_{SPH}\right)$')
        self.line, = self.ax.plot(self.t,
                                  self.p,
                                  label=r'$p_{SPH}$',
                                  color="black",
                                  linewidth=1.0)
        # Set some options
        self.ax.grid()
        self.ax.legend(loc='best')
        self.ax.set_autoscale_on(False)
        self.ax.set_xlabel(r"$t \, [\mathrm{s}]$", fontsize=21)
        self.ax.set_ylabel(r"$p \, [\mathrm{Pa}]$", fontsize=21)
        self.ax.set_xlim(0, 5)
        self.ax.set_ylim(-1000, 5000)
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
        self.t = data[0]
        self.p = []
        self.pmin = []
        self.pmax = []
        for i in range(len(data[0])):
            p = data[3][i]
            rho = data[4][i]
            rho_var = data[5][i]
            pmin = eos(rho - rho_var)
            pmax = eos(rho + rho_var)
            sum_w = data[6][i]
            # Two filters may help to get better pressure records:
            #  * Small kernels completions, which means not enough particles
            #    to get a confiable point
            #  * Large data deviation, which usually means that only free
            #    surface particles are computed.
            if sum_w < 0.01:
                self.p.append(0.0)
                self.pmin.append(0.0)
                self.pmax.append(0.0)
                continue
            if rho_var / refd > 0.002:
                self.p.append(0.0)
                self.pmin.append(0.0)
                self.pmax.append(0.0)
                continue
            self.p.append(p)
            self.pmin.append(pmin)
            self.pmax.append(pmax)
            # self.rho_var.append(sum_w)
        for coll in (self.ax.collections):
            self.ax.collections.remove(coll)


        self.line.set_data(self.t, self.p)
        # print(dir(self.var_line))
        # print(self.var_line.get_paths())
        # self.var_line.set_data(self.t, self.pmin, self.pmax)
        # Redraw
        self.var_line = self.ax.fill_between(self.t,
                                             self.pmin,
                                             self.pmax,
                                             interpolate=True,
                                             facecolor=(0.75, 0.75, 0.75, 1.0),
                                             edgecolor=(0.5, 0.5, 0.5, 0.0),
                                             linestyle='--',
                                             label=r'$p_{SPH} \pm \sigma\left(p_{SPH}\right)$')
        fake_artist = plt.Rectangle((0, 0), 1, 1,
                                    facecolor=(0.75, 0.75, 0.75, 1.0),
                                    edgecolor=(1.0, 1.0, 1.0, 0.0))
        self.ax.legend([self.exp_line,
                        self.line,
                        fake_artist],
                       [self.exp_line.get_label(),
                        self.line.get_label(),
                        self.var_line.get_label()],
                       loc='best')
        # self.ax.legend(loc='best')
        self.fig.canvas.draw()


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    widget = FigureController()
    widget.setWindowTitle("Sensor 1")
    widget.show()
    sys.exit(app.exec_())
