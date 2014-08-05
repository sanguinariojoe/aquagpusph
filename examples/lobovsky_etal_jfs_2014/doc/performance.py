#########################################################
#                                                       #
#    #    ##   #  #   #                           #     #
#   # #  #  #  #  #  # #                          #     #
#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #
#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #
#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #
#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #
#                            # #             #          #
#                          ##  #             #          #
#                                                       #
#########################################################
#
# timestep.py
#
# Python script to read the timesteps and the time
# consumed to plot the device performance.
#
#
#########################################################

import sys
import os
try:
	from PyQt4 import QtGui
except:
	print("PyQt4 is required to can use this tool.")
	sys.exit(255)

try:
	from matplotlib.figure import Figure
	from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
except:
	print("matplotlib is required to can use this tool.")
	sys.exit(255)

# Log file path
log_id = 0
log_path = 'log.{}.html'.format(log_id)
while os.path.isfile('log.{}.html'.format(log_id + 1)):
    log_id += 1
    log_path = 'log.{}.html'.format(log_id)

class FigureController(FigureCanvas):
	"""Matplotlib figure widget controller"""
	def __init__(self):
		"""Constructor"""
		# Create the figure in the canvas
		self.fig = Figure()
		self.ax = self.fig.add_subplot(111)
		FigureCanvas.__init__(self, self.fig)
		# generates first "empty" plot
		self.performance = [0.0]
		self.time = [0.0]
		self.line, = self.ax.plot(self.time,
		                          self.performance,
		                          label='Performance',
		                          color="black",
		                          linewidth=2.0)
		# Set some options
		self.ax.grid()
		self.ax.set_autoscale_on(False)
		self.ax.set_xlabel(r"$t_{sim} \, [\mathrm{s}]$", fontsize=21)
		self.ax.set_ylabel(r"$\frac{t_{CPU}}{step} \, [\mathrm{s}]$", fontsize=21)
		# force the figure redraw
		self.fig.canvas.draw()
		# call the update method (to speed-up visualization)
		self.timerEvent(None)
		# start timer, trigger event every 1000 millisecs (=1sec)
		self.timer = self.startTimer(1000)

	def readLines(self, date_line, sim_line):
		""" Extract the data from the lines """
		# Extract the date
		index = date_line.find("<i>")+3
		line = date_line[index:]
		index = line.find("Printed file")-2
		d = line[:index]
		# Extract the step
		index = sim_line.find("nstep=")+6
		line = sim_line[index:]
		index = line.find("n=")-2
		s = line[:index]
		# Extract the time
		index = sim_line.find("time=")+5
		line = sim_line[index:]
		index = line.find("dt=")-2
		t = line[:index]
		return (d,s,t)
		
	def readDate(self, date):
		""" Extract the date day, hour, minute, second """
		# Clear double separators
		edit = date[:]
		while(edit.find("  ") != -1):
			edit = edit.replace("  ", " ")
		# Get the date fields:
		#   -# Day on week
		#   -# Month (by name)
		#   -# Day
		#   -# HH:MM:SS
		#   -# Year
		fields = edit.split(" ")
		day = int(fields[2])
		# Get the hour fields:
		#   -# Hour
		#   -# Minute
		#   -# Second
		fields = fields[3].split(":")
		hour = int(fields[0])
		minute = int(fields[1])
		second = int(fields[2])
		return (day,hour,minute,second)

	def readFile(self):
		""" Read and extract data from the log file """
		# Read the file by lines
		f = open(log_path, "r")
		lines = f.readlines()
		f.close()
		# Run on the lines looking for "Printed file" sentences
		dates = []
		steps = []
		times = []
		for i,l in enumerate(lines):
			if "Printed file" in l:
				# The lines must be analised
				d,s,t = self.readLines(lines[i], lines[i+1])
				dates.append(d)
				steps.append(s)
				times.append(t)
		# Decript the data
		self.time = []
		self.performance = []
		for i in range(1,len(dates)):
			day1,hour1,minute1,second1 = self.readDate(dates[i])
			day0,hour0,minute0,second0 = self.readDate(dates[i-1])
			seconds  = (day1-day0)*24*3600
			seconds += (hour1-hour0)*3600
			seconds += (minute1-minute0)*60
			seconds += second1-second0
			self.time.append(0.5*(float(times[i]) + float(times[i-1])))
			self.performance.append(seconds/(float(steps[i])-float(steps[i-1])))

	def timerEvent(self, evt):
		"""Custom timerEvent code, called at timer event receive"""
		# Read the log file
		self.readFile()
		if len(self.time) <= 1:
			return
		# Set the new data and limits
		self.line.set_data(self.time, self.performance)
		self.ax.set_xlim(0, self.time[-1])
		min_val = min(self.performance)
		max_val = max(self.performance)
		self.ax.set_ylim(min_val - 0.1*abs(min_val),
		                 max_val + 0.1*abs(max_val))
		# Redraw
		self.fig.canvas.draw()

# create the Application
app = QtGui.QApplication(sys.argv)
widget = FigureController()
widget.setWindowTitle("Performance along the simulation")
widget.show()
# And start working
sys.exit(app.exec_())

