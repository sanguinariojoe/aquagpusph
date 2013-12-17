#!/bin/bash

if [[ $1 == "--run" ]]; then
	@BINARY_DIR@/AQUAgpusph -i @EXAMPLE_DEST_DIR@/Main.xml --no-reassembly
elif [[ $1 == "--plot" ]]; then
	if [ ! -f zMovement.dat ]; then
		echo ""
		echo "zMovement.dat output file is not yet ready!"
		echo ""
		echo "    You must execute run.sh --plot on the same folder where run.sh --run has been launched"
		echo ""
		exit 255
	fi
	gnuplot @EXAMPLE_DEST_DIR@/doc/plot.gnuplot
else
	echo ""
	echo "Usage: run.sh [--run/--plot]"
	echo ""
	echo "run.sh --run"
	echo "    Runs the simulation. You can stop simulation pressing 'c' key"
	echo "run.sh --plot"
	echo "    Plots in real time the output results"
	echo ""
fi
