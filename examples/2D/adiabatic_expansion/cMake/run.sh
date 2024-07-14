#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Fluid.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f output.*.vtu
    rm -f output.vtu.series
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph -d 2 -i Main.xml
elif [[ $1 == "--plot" ]]; then
    python @EXAMPLE_DEST_DIR@/plot_$2.py
else
    echo ""
    echo "Usage: run.sh [--run/--plot x/--plot e/--plot f]"
    echo ""
    echo "run.sh --run"
    echo "    Runs the simulation"
    echo "run.sh --plot x"
    echo "    Plot the piston motion"
    echo "run.sh --plot e"
    echo "    Plot the resulting energy"
    echo "run.sh --plot f"
    echo "    Plot the resulting force on the wall"
    echo ""
fi
