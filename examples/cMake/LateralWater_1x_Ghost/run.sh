#!/bin/bash

if [[ $1 == "--run" ]]; then
	rm -f Sensors.dat
	@BINARY_DIR@/AQUAgpusph2D -i @EXAMPLE_DEST_DIR@/Main.xml --no-reassembly
elif [[ $1 == "--plot" ]]; then
	if [ ! -f Sensors.dat ]; then
		echo ""
		echo "Failure while opening Sensors.dat output file"
		echo ""
		echo " * You must execute run.sh --plot on the same folder where run.sh --run has been launched"
		echo " * AQUAgpusph may take some time before generate this file"
		echo ""
		exit 255
	fi
	cp -f @EXAMPLE_DEST_DIR@/doc/lateral_water_1x.txt lateral_water_1x.txt
	gnuplot @EXAMPLE_DEST_DIR@/doc/plot.gnuplot
else
	echo ""
	echo "Usage: run.sh [--run/--plot]"
	echo ""
	echo "run.sh --run"
	echo "    Runs the simulation. You can stop simulation pressing the 'c' key"
	echo "run.sh --plot"
	echo "    Plots in real time the output results"
	echo ""
fi
