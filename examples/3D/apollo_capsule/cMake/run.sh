#!/bin/bash

set -e

if [[ $1 == "--run" ]]; then
    rm -f Titan.dat
    rm -f Pool.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f apollo*.vtu
    rm -f apollo.vtu.series
    rm -f pool*.vtu
    rm -f pool.vtu.series
    rm -f Timing.dat
    rm -f force*.out
    rm -f motion*.out
    rm -f Performance.dat
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph -i Main.xml
elif [[ $1 == "--plot" ]]; then
    python @EXAMPLE_DEST_DIR@/plot_$2.py
else
    echo ""
    echo "Usage: run.sh [--run/--plot]"
    echo ""
    echo "run.sh --plot m"
    echo "    Plot the resulting motion"
    echo "run.sh --plot t"
    echo "    Plot the performance"
    echo ""
fi
