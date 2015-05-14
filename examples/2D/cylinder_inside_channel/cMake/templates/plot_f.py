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


# Input data from the case generation script
rho = {{REFD}}
D = {{D}}
U = {{U}}
# Coefficients adimensionalization factor
COEFF_FAC = 1.0 / (0.5 * rho * D * U**2)
TIME_FAC = {{U}} / {{D}}


class FigureController(FigureCanvas):
    """Matplotlib figure widget controller"""

    def __init__(self):
        """Constructor"""
        # Create the figure in the canvas
        self.fig = Figure()
        self.ax_cl = self.fig.add_subplot(111)
        self.ax_cd = self.ax_cl.twinx()
        FigureCanvas.__init__(self, self.fig)
        # generates first "empty" plot
        self.cl, = self.ax_cl.plot([0.0],
                                   [0.0],
                                   label=r'$c_L$',
                                   color="black",
                                   linestyle="-",
                                   linewidth=1.0)
        self.cd, = self.ax_cd.plot([0.0],
                                   [0.0],
                                   label=r'$c_D$',
                                   color="red",
                                   linestyle="-",
                                   linewidth=1.0)
        # Set some options
        self.ax_cl.grid()
        self.ax_cl.set_xlim(0, 1.0)
        self.ax_cl.set_ylim(-1.0, 1.0)
        self.ax_cl.set_autoscale_on(False)
        self.ax_cl.set_xlabel(r"$t \, U / D$", fontsize=21)
        self.ax_cl.set_ylabel(r"$c_L$", fontsize=21)
        self.ax_cd.grid()
        self.ax_cd.set_xlim(0, 1.0)
        self.ax_cd.set_ylim(-1.0, 1.0)
        self.ax_cd.set_autoscale_on(False)
        self.ax_cd.set_ylabel(r"$c_D$", fontsize=21, color="red")
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
        data = self.readFile('PressureForces.dat')
        t = data[0]
        fpx = data[1]
        fpy = data[2]
        data = self.readFile('ViscousForces.dat')
        t = data[0]
        fvx = data[1]
        fvy = data[2]
        if(len(tv) < len(tp)):
            tp = tp[0:len(tv)]
            fpx = fpx[0:len(tv)]
            fpy = fpy[0:len(tv)]
        elif(len(tv) > len(tp)):
            tv = tv[0:len(tp)]
            fvx = fvx[0:len(tp)]
            fvy = fvy[0:len(tp)]
        t = []
        fx = []
        fy = []
        for i in range(len(tp)):
            t.append(tp[i] * TIME_FAC)
            fx.append((fpx[i] + fvx[i]) * COEFF_FAC)
            fy.append((fpy[i] + fvy[i]) * COEFF_FAC)
        self.cl.set_data(t, fy)
        self.cd.set_data(t, fx)

        clmax = max(max(fy), abs(min(fy)))
        cdmax = max(max(fx), abs(min(fx)))
        self.ax_cl.set_xlim(0, 1.05 * max(t))
        self.ax_cd.set_xlim(0, 1.05 * max(t))
        self.ax_cl.set_ylim(-1.05 * clmax, 1.05 * clmax)
        self.ax_cd.set_ylim(-1.05 * cdmax, 1.05 * cdmax)
        for tl in self.ax_cd.get_yticklabels():
            tl.set_color("red")

        # Redraw
        self.fig.canvas.draw()


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    widget = FigureController()
    widget.setWindowTitle("Lift and drag coefficients")
    widget.show()
    sys.exit(app.exec_())
