#!/bin/bash
# Submit all neurodamus jobs in refitting_results/fitting/n100
# Called as a SLURM job after inject_thresholds.py completes
# Usage: sbatch --dependency=afterok:INJECT_JOB submit_after_inject.sh [freq_filter]

REFDIR="/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100"
FREQ_FILTER="${1:-}"  # optional, e.g. "10Hz"

submitted=0
skipped=0

find "$REFDIR" -name "neurodamus_sbatch.sh" | sort | while read sbatch_file; do
  dir=$(dirname "$sbatch_file")
  sim_label=$(basename "$dir")
  
  # Apply frequency filter if set
  if [ -n "$FREQ_FILTER" ] && [[ "$sim_label" != *"$FREQ_FILTER"* ]]; then
    continue
  fi
  
  # Skip if already has a SUCCESS marker
  if ls "$dir"/simulation_config.json.SUCCESS &>/dev/null; then
    echo "SKIP (done): $dir"
    continue
  fi
  
  pushd "$dir" > /dev/null
  jid=$(sbatch --parsable "$sbatch_file")
  echo "SUBMITTED $jid: $dir"
  popd > /dev/null
done

echo "All done."
