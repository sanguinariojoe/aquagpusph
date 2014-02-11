#!/bin/bash

function sedeasy {
	sed -i "s/$(echo $1 | sed -e 's/\([[\/.*]\|\]\)/\\&/g')/$(echo $2 | sed -e 's/[\/&]/\\&/g')/g" $3
}

if [[ $1 == "--run" ]]; then
	rm -f Sensors.dat
	@BINARY_DIR@/AQUAgpusph2D -i @EXAMPLE_DEST_DIR@/Main.xml --no-reassembly
elif [[ $1 == "--plot" ]]; then
	if [ ! -f Sensors.dat ]; then
		echo ""
		echo "Sensors.dat output file is not ready yet!"
		echo ""
		echo "    You must execute run.sh --plot on the same folder where run.sh --run has been launched"
		echo ""
		exit 255
	fi
	# Replace the path doc/ by the absolute files path
	sedeasy "@EXAMPLE_DEST_DIR@/doc/" "doc/" @EXAMPLE_DEST_DIR@/doc/$2.gnuplot
	sedeasy "doc/" "@EXAMPLE_DEST_DIR@/doc/" @EXAMPLE_DEST_DIR@/doc/$2.gnuplot
	gnuplot @EXAMPLE_DEST_DIR@/doc/$2.gnuplot
elif [[ $1 == "--performance" ]]; then
	if [ ! -f Log.html ]; then
		echo ""
		echo "Log.html file is not ready yet!"
		echo ""
		echo "    You must execute run.sh --timestep on the same folder where run.sh --run has been launched"
		echo ""
		exit 255
	fi
	python @EXAMPLE_DEST_DIR@/doc/performance.py
else
	echo ""
	echo "Usage: run.sh [--run/--plot/--timestep] [options...]"
	echo ""
	echo "run.sh --run"
	echo "    Runs the simulation. You can stop simulation pressing 'c' key"
	echo "run.sh --performance"
	echo "    Plots in real time the performance (timesteps per second)"
	echo "run.sh --plot [graph]"
	echo "    Plots in real time the output results. Valid graph are:"
	echo "    sensor1    1st pessure sensor"
	echo "    sensor2    2nd pessure sensor"
	echo "    sensor3    3rd pessure sensor"
	echo "    sensor4    4th pessure sensor"
	echo "    sensors    All the pessure sensors simultaneously"
	echo "    x          Propagation of the surge front"
	echo ""
fi
