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
        self.ax11 = self.fig.add_subplot(221)
        self.ax21 = self.fig.add_subplot(222)
        self.ax12 = self.fig.add_subplot(223)
        self.ax22 = self.fig.add_subplot(224)
        self.ax = (self.ax11, self.ax21, self.ax12, self.ax22)
        FigureCanvas.__init__(self, self.fig)

        FNAME = path.join('@EXAMPLE_DEST_DIR@', 'test_case_2_exp_data.dat')
        # For some reason the input file is bad sortened
        T,_,_,_,_,_,_,_,_,H3,H2,H1,H4, = self.readFile(FNAME)
        exp_t = T
        exp_h = (H1, H2, H3, H4)
        titles = ('H1', 'H2', 'H3', 'H4')
        self.lines = []
        for i in range(len(self.ax)):
            ax = self.ax[i]
            t = [0.0]
            h = [0.0]
            line, = ax.plot(t,
                            h,
                            label=r'$H_{SPH}$',
                            color="black",
                            linewidth=1.0)
            self.lines.append(line)
            ax.plot(exp_t,
                    exp_h[i],
                    label=r'$H_{Exp}$',
                    color="red",
                    linewidth=1.0)
            # Set some options
            ax.grid()
            ax.legend(loc='best')
            ax.set_title(titles[i])
            ax.set_xlim(0, 6)
            ax.set_ylim(0.0, 0.6)
            ax.set_autoscale_on(False)
            ax.set_xlabel(r"$t \, [\mathrm{s}]$", fontsize=21)
            ax.set_ylabel(r"$H \, [\mathrm{m}]$", fontsize=21)
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
        data = self.readFile('sensors_h.out')
        t = data[0]
        hh = (data[-4], data[-3], data[-2], data[-1])
        for i in range(len(hh)):
            h = hh[i]
            self.lines[i].set_data(t, h)

        # Redraw
        self.fig.canvas.draw()


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    widget = FigureController()
    widget.setWindowTitle("Wave height")
    widget.show()
    sys.exit(app.exec_())
