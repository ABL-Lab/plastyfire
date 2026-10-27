#!/bin/bash
source /home/dhuruva/ws/plastyfire/setupenv.sh
python /home/dhuruva/ws/plastyfire/submit_stdp_traces.py --max-pairs 2 --execution-mode cpu --workers 30
