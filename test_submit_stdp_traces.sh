#!/bin/bash
source /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/setupenv.sh
python /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/submit_stdp_traces.py --max-pairs 2 --execution-mode cpu --workers 30
