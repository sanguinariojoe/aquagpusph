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
else
    echo ""
    echo "Usage: run.sh --run"
    echo ""
fi
