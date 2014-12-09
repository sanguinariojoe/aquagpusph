#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Sensors.dat
    rm -f Fluid.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f output.*.vtu
    rm -f output.pvd
    rm -f sensors.*.vtu
    rm -f sensors.pvd
    rm -f sensors.out
    rm -f Timing.dat
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph -i Main.xml
elif [[ $1 == "--plot" ]]; then
    if [ ! -f sensors.out ]; then
        echo ""
        echo "Failure while opening sensors.out output file"
        echo ""
        echo "run.sh --plot should be executed on the same folder of run.sh --run"
        echo ""
        exit 255
    fi
    python @EXAMPLE_DEST_DIR@/plot_p.py
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
