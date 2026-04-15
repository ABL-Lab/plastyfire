#!/bin/bash
# Submit all neurodamus_sbatch.sh scripts under refitting_results.
# Usage:
#   ./submit_neurodamus.sh                  # submit all frequencies
#   ./submit_neurodamus.sh --freq 10        # only 10Hz
#   ./submit_neurodamus.sh --freq 10 --dry-run

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/refitting_results"

FREQ_FILTER=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --freq)   FREQ_FILTER="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        *) echo "Unknown option: $1"; echo "Usage: $0 [--freq HZ] [--dry-run]"; exit 1 ;;
    esac
done

if [[ -n "$FREQ_FILTER" ]]; then
    PATTERN="${FREQ_FILTER}Hz_*"
else
    PATTERN="*Hz_*"
fi

mapfile -t SCRIPTS < <(find "$RESULTS_DIR" -path "*/$PATTERN/neurodamus_sbatch.sh" | sort)

if [[ ${#SCRIPTS[@]} -eq 0 ]]; then
    echo "No neurodamus_sbatch.sh scripts found matching freq='${FREQ_FILTER:-any}' under $RESULTS_DIR"
    exit 0
fi

echo "Found ${#SCRIPTS[@]} script(s) [freq filter: ${FREQ_FILTER:-any}]"
[[ $DRY_RUN -eq 1 ]] && echo "(dry-run: not submitting)"

SUBMITTED=0
for SCRIPT in "${SCRIPTS[@]}"; do
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "  [dry-run] sbatch $SCRIPT"
    else
        JOB=$(sbatch --chdir="$(dirname "$SCRIPT")" "$SCRIPT" 2>&1)
        echo "  $JOB  <- $SCRIPT"
        ((SUBMITTED++))
    fi
done

[[ $DRY_RUN -eq 0 ]] && echo "Submitted $SUBMITTED job(s)."
