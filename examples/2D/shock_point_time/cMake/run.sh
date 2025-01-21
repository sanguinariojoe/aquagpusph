#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Fluid.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f output.*.vtu
    rm -f *.vtu
    rm -f output.vtu.series
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph -d 2 -i Main.xml
else
    echo ""
    echo "Usage: run.sh [--run/--plot e/--plot f]"
    echo ""
    echo "run.sh --run"
    echo "    Runs the simulation. You can stop simulation pressing the 'c' key"
    echo "run.sh --plot e"
    echo "    Plot the resulting energy"
    echo "run.sh --plot f"
    echo "    Plot the resulting force on the wall"
    echo "run.sh --plot r"
    echo "    Plot the residues statistics"
    echo ""
fi