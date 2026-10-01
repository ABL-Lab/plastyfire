#!/bin/bash
#SBATCH --job-name=offload_ljp25
#SBATCH --account=rrg-emuller
# Inode quota relief: move the REJECTED ljp25 raw prefire outputs (bluecellulab_results_delta-ljp25-prefire-vseg, 3064 files,
# 948 already moved) to /scratch/dhuruva/plastyfire_offload (same relative paths; move back with rsync to undo).
# Extracted ljp25 npz stay on project. MEASURED 22127439: TIMEOUT at 0:30, 1.99 GB MaxRSS (of 2G), 83% CPU; rsync of large pkl over lustre is slow -> 3G, 2:00.
#   Record seff here after the run.
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=02:00:00
#SBATCH --chdir=/project/rrg-emuller/dhuruva/plastyfire
#SBATCH --output=/scratch/dhuruva/plastyfire_offload/offload_ljp25_%j.out
set -uo pipefail
O=/scratch/dhuruva/plastyfire_offload
find refitting_results -path '*bluecellulab_results_delta-ljp25-prefire-vseg/*' -type f > $O/ljp25_left.txt
echo "to move: $(wc -l < $O/ljp25_left.txt)"
rsync -a --files-from=$O/ljp25_left.txt --remove-source-files . $O/
find refitting_results -type d -name 'bluecellulab_results_delta-ljp25-prefire-vseg' -empty -delete
left=$(find refitting_results -path '*bluecellulab_results_delta-ljp25-prefire-vseg/*' -type f | wc -l)
echo "left on project: $left; on scratch: $(find $O/refitting_results -type f | wc -l)"
[ "$left" -eq 0 ]
