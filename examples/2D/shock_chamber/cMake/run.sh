#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Fluid.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f output.*.vtu
    rm -f output.vtu.series
    rm -f Timing.dat
    rm -f Performance.dat
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph -d 2 -i Main.xml
elif [[ $1 == "--plot" ]]; then
    python @EXAMPLE_DEST_DIR@/plot_$2.py
else
    echo ""
    echo "Usage: run.sh [--run/--plot p/--plot h]"
    echo ""
    echo "run.sh --run"
    echo "    Runs the simulation. You can stop simulation pressing the 'c' key"
    echo "run.sh --plot t"
    echo "    Plot the performance"
    echo ""
fi
