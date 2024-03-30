#!/bin/bash

if [[ $1 == "--run" ]]; then
    rm -f Fluid.dat
    rm -f BC.dat
    rm -f AQUAgpusph.save.*.xml
    rm -f log.*.html
    rm -f *.vtu
    rm -f *.pvd
    rm -f Timing.dat
    rm -f Performance.dat
    @EXAMPLE_DEST_DIR@/Create.py
    @BINARY_DIR@/AQUAgpusph2D -i Main.xml
else
    echo ""
    echo "Usage: run.sh --run"
    echo ""
fi
