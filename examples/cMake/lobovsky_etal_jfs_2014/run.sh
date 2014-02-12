#!/bin/bash

if [[ $1 == "--run" ]]; then
	rm -f Sensors.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f energy.*.dat
    rm -f log.*.html
    rm -f output.*.vtu
    rm -f output.pvd
	@BINARY_DIR@/AQUAgpusph2D -i @EXAMPLE_DEST_DIR@/Main.xml
elif [[ $1 == "--continue" ]]; then
    INPUT=`ls -1 @EXAMPLE_DEST_DIR@ | grep AQUAgpusph.save. | tail -n 1`
    if [[ $INPUT == "" ]]; then
		echo ""
		echo "A backup file AQUAgpusph.save.*.xml cannot be found"
		echo ""
		exit 255
	fi
    @BINARY_DIR@/AQUAgpusph2D -i $INPUT
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
	python @EXAMPLE_DEST_DIR@/doc/$2.py
elif [[ $1 == "--performance" ]]; then
	if [ ! -f Log.html ]; then
		echo ""
		echo "Failure while opening Log.html output file"
		echo ""
		echo " * You must execute run.sh --performance on the same folder where run.sh --run has been launched"
		echo " * AQUAgpusph may take some time before generate this file"
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
	echo ""
fi
