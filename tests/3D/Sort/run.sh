#!/bin/bash

set -e

# Run the "simulation"
../../../bin/AQUAgpusph -i main.xml

# Check the result
python check.py

# Clear everything
rm -f AQUAgpusph.save.* log.* set0.*.dat set1.*.dat
