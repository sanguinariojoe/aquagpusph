#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Sensors.dat
    rm -f Fluid.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f output.*.vtu
    rm -f output.vtu.series
    rm -f sensors.*.vtu
    rm -f sensors.vtu.series
    rm -f sensors.out
    rm -f sensors_h.out
    rm -f Timing.dat
    rm -f Performance.dat
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/@BINARY_NAME@ -i Main.xml
elif [[ $1 == "--plot" ]]; then
    python @EXAMPLE_DEST_DIR@/plot_$2.py
else
    echo ""
    echo "Usage: run.sh [--run/--plot]"
    echo ""
    echo "run.sh --run"
    echo "    Runs the simulation. You can stop simulation pressing the 'c' key"
    echo "run.sh --plot p"
    echo "    Plot the pressure"
    echo "run.sh --plot h"
    echo "    Plot the wave height"
    echo "run.sh --plot r"
    echo "    Plot the residues statistics"
    echo ""
fi
