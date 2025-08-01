#!/bin/bash

set -e

rm -f AQUAgpusph.save.* log.* output.*.dat *.data performance.json
@BINARY_DIR@/@BINARY_NAME@ -d 2 -i main.xml

# Check the result
python @TEST_DEST_DIR@/check.py
