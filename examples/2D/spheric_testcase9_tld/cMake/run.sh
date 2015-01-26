#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Fluid.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f output.*.vtu
    rm -f output.pvd
    rm -f Timing.dat
    rm -f Forces.dat
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph2D -i Main.xml
elif [[ $1 == "--plot" ]]; then
    if [ ! -f Motion.dat ]; then
        echo ""
        echo "Failure while opening Motion.dat output file"
        echo ""
        echo " * You must execute run.sh --plot on the same folder where run.sh --run has been launched"
        echo " * AQUAgpusph may take some time before generate this file"
        echo ""
        exit 255
    fi
    python @EXAMPLE_DEST_DIR@/plot.py
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