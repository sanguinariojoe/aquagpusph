#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Fluid.dat
    rm -f Tank.dat
    rm -f FreeSurface.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f Fluid.*.vtu
    rm -f Fluid.pvd
    rm -f FreeSurface.*.vtu
    rm -f FreeSurface.pvd
    rm -f output.*.vtu
    rm -f output.pvd
    rm -f output_bc.*.vtu
    rm -f output_bc.pvd
    rm -f Timing.dat
    rm -f Forces.dat
    rm -f Performance.dat
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph2D -i particlesPacking.xml
    rm -f FreeSurface.*.vtu
    rm -f FreeSurface.pvd
    @BINARY_DIR@/AQUAgpusph2D -i Main.xml
    rm -f Fluid.*.vtu
    rm -f Fluid.pvd
elif [[ $1 == "--plot" ]]; then
    python @EXAMPLE_DEST_DIR@/plot_$2.py
else
    echo ""
    echo "Usage: run.sh [--run/--plot m/--plot t]"
    echo ""
    echo "run.sh --run"
    echo "    Runs the simulation. You can stop simulation pressing the 'c' key"
    echo "run.sh --plot m"
    echo "    Plot the resulting motion"
    echo "run.sh --plot t"
    echo "    Plot the performance"
    echo ""
fi
