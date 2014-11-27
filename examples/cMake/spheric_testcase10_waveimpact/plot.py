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
        self.ax = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
        # generates first "empty" plot
        FNAME = path.join('@EXAMPLE_DEST_DIR@', 'lateral_water_1x.txt')
        T,P,A,DADT,_,_ = np.loadtxt(FNAME,
                                    delimiter='\t',
                                    skiprows=1,
                                    unpack=True)
        self.exp_t = T
        self.exp_p = 100.0 * P
        self.exp_line, = self.ax.plot(self.exp_t,
                                      self.exp_p,
                                      label=r'$p_{Exp}$',
                                      color="red",
                                      linewidth=1.0)
        self.t = [0.0]
        self.p = [0.0]
        self.line, = self.ax.plot(self.t,
                                  self.p,
                                  label=r'$p_{SPH}$',
                                  color="black",
                                  linewidth=1.0)
        # Set some options
        self.ax.grid()
        self.ax.legend(loc='best')
        self.ax.set_xlim(0, 5)
        self.ax.set_ylim(-1000, 5000)
        self.ax.set_autoscale_on(False)
        self.ax.set_xlabel(r"$t \, [\mathrm{s}]$", fontsize=21)
        self.ax.set_ylabel(r"$p \, [\mathrm{Pa}]$", fontsize=21)
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
        for l in lines[1:]:
            l = l.strip()
            while l.find('  ') != -1:
                l = l.replace('  ', ' ')
            fields = l.split(' ')
            try:
                data.append(map(float, fields))
            except:
                continue
        # Transpose the data
        return map(list, zip(*data))

    def timerEvent(self, evt):
        """Custom timerEvent code, called at timer event receive"""
        # Read and plot the new data
        data = self.readFile('sensors.out')
        self.t = data[0]
        self.p = data[1]
        self.line.set_data(self.t, self.p)

        # Redraw
        self.fig.canvas.draw()


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    widget = FigureController()
    widget.setWindowTitle("Pressure")
    widget.show()
    sys.exit(app.exec_())
