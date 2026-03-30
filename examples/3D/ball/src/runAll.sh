#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python3 $SCRIPT_DIR/ball_new.py
chmod +x clean.sh
chmod +x run.sh
./clean.sh
./run.sh
