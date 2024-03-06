#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Fluid.dat
    rm -f Sensors.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f output.*.vtu
    rm -f sensors.*.vtu
    rm -f output.pvd
    rm -f sensors.pvd
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph2D -i Main.xml
elif [[ $1 == "--plot" ]]; then
    python @EXAMPLE_DEST_DIR@/plot_$2.py
else
    echo ""
    echo "Usage: run.sh [--run/--plot e]"
    echo ""
    echo "run.sh --run"
    echo "    Runs the simulation"
    echo "run.sh --plot e"
    echo "    Plot the resulting energy"
    echo ""
fi
