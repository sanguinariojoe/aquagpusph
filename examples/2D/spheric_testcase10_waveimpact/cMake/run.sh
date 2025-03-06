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
    rm -f Timing.dat
    rm -f force.out
    rm -f Performance.dat
    rm -f Energy.dat
    rm -f midpoint.out
    python Create.py
    @BINARY_DIR@/@BINARY_NAME@ -d 2 -i Main.xml
elif [[ $1 == "--plot" ]]; then
    python @EXAMPLE_DEST_DIR@/plot_$2.py
else
    echo ""
    echo "Usage: run.sh [--run/--plot [PLOTTER]]"
    echo ""
    echo "run.sh --run"
    echo "    Run the simulation"
    echo "run.sh --plot p"
    echo "    Plot the pressure in the sensor"
    echo "run.sh --plot e"
    echo "    Plot the energy profile"
    echo "run.sh --plot r"
    echo "    Plot the residues statistics"
    echo ""
fi
