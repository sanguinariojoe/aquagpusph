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
import sys
import os
from os import path
try:
    from PyQt4 import QtGui, QtCore
except:
    try:
        from PySide import QtGui, QtCore
    except:
        raise ImportError("PyQt4 or PySide is required to use this tool")

try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
except:
    raise ImportError("matplotlib is required to use this tool")


class FigureController(FigureCanvas):
    """Matplotlib figure widget controller"""

    def __init__(self, parent=None):
        """Constructor"""
        # Create the figure in the canvas
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        # generates first "empty" plot
        t = [0.0]
        e = [0.0]
        self.fove = self.ax.fill_between(t,
                                         0,
                                         e,
                                         facecolor='red',
                                         linewidth=0.0)
        self.fave = self.ax.fill_between(t,
                                         0,
                                         e,
                                         facecolor='blue',
                                         linestyle="-",
                                         linewidth=0.0)
        self.love, = self.ax.plot(t,
                                  e,
                                  color='#990000',
                                  linestyle="-",
                                  linewidth=2.0,
                                  label='Average overhead')
        self.lave, = self.ax.plot(t,
                                  e,
                                  color='#000099',
                                  linestyle="-",
                                  linewidth=2.0,
                                  label='Average elapsed')
        self.line, = self.ax.plot(t,
                                  e,
                                  color="black",
                                  linestyle="-",
                                  linewidth=1.0,
                                  alpha=0.5,
                                  label='Elapsed')
        # Set some options
        self.ax.grid()
        self.ax.set_xlim(0, 0.1)
        self.ax.set_ylim(-0.1, 0.1)
        self.ax.set_autoscale_on(False)
        self.ax.set_xlabel(r"$t \, [\mathrm{s}]$", fontsize=21)
        self.ax.set_ylabel(r"$t_{CPU} \, [\mathrm{s}]$", fontsize=21)
        self.ax.legend(handles=[self.lave, self.love, self.line],
                       loc='upper right')
        # force the figure redraw
        self.fig.canvas.draw()
        # call the update method (to speed-up visualization)
        self.timerEvent(None)
        # start timer, trigger event every 10000 millisecs (=10sec)
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
        data = self.readFile('Performance.dat')
        t = data[0]
        e = data[1]
        e_ela = data[2]
        e_ove = data[5]
        # Clear nan values
        for i in range(len(e_ela)):
            if math.isnan(e_ela[i]):
                e_ela[i] = 0.0
            if math.isnan(e_ove[i]):
                e_ove[i] = 0.0
        e_ave = [e_ela[i] - e_ove[i] for i in range(len(e_ela))]
        # clear the fills
        for coll in (self.ax.collections):
            self.ax.collections.remove(coll)
        self.fove = self.ax.fill_between(t,
                                         0,
                                         e_ela,
                                         facecolor='red',
                                         linestyle="-",
                                         linewidth=2.0)
        self.fave = self.ax.fill_between(t,
                                         0,
                                         e_ave,
                                         facecolor='blue',
                                         linestyle="-",
                                         linewidth=2.0)
        self.love.set_data(t, e_ela)
        self.lave.set_data(t, e_ave)
        self.line.set_data(t, e)

        self.ax.set_xlim(0, t[-1])
        self.ax.set_ylim(0, 1.5 * e_ela[-1])

        # Redraw
        self.fig.canvas.draw()


class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        self.main_widget = QtGui.QWidget(self)

        lv = QtGui.QVBoxLayout(self.main_widget)
        plot_widget = FigureController(self.main_widget)
        show_inst_ela = QtGui.QCheckBox('Show instantaneous elapsed time')
        show_inst_ela.setCheckState(True)
        show_inst_ela.setTristate(False)
        lv.addWidget(plot_widget)
        lv.addWidget(show_inst_ela)

        show_inst_ela.stateChanged.connect(self.showInstEla)

        self.show_inst_ela = show_inst_ela
        self.plot_widget = plot_widget

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

    def showInstEla(self, state):
        if not state:
            self.plot_widget.line.set_alpha(0.0)
            self.plot_widget.ax.legend(handles=[self.plot_widget.lave,
                                                self.plot_widget.love],
                                       loc='upper right')

            return
        self.plot_widget.ax.legend(handles=[self.plot_widget.lave,
                                            self.plot_widget.love,
                                            self.plot_widget.line],
                                   loc='upper right')
        self.plot_widget.line.set_alpha(0.5)
        return

    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    aw = ApplicationWindow()
    aw.setWindowTitle("Performance")
    aw.show()
    sys.exit(app.exec_())