#!/bin/bash
# Representative traces of fit v5_S1_V5c_s6_39 (v5c, emodel delta-split1): BCL live (live_v5.run_task + recording hook)
# vs the offline kernel, 11 records (traces_v5.GROUPS) -> /scratch/dhuruva/bcl_traces_S1_V5c/*.npz; then the figures.
# STAGE=traces (default): traces_v5.py --check, then the run (11 workers + main).
# STAGE=plot: plot_traces_v5.py -> results/traces_S1_V5c_<group>.png/.pdf + overview; submit with
#   sbatch --cpus-per-task=1 --mem=4G --time=00:15:00 (afterok the traces job).
# Sizing basis: BCL pilot 22144694 MEASURED 8:54 for 38 records at 12 workers, 5.24 GB (~3 core-min, ~0.4 GB per
#   record); 11 records -> 12 CPU, 8G. MEASURED 22148298: 5:50, 5.88 GB (np.array(Vector) crash, 0 npz) -> 0:15.
#   MEASURED 22148840: 6:18, 7.28 GB/8G, CPU 44% at 12 (records 92-329 s, ~34 core-min) -> 7 CPU (6 workers), 8G, 0:15.
#SBATCH --job-name=traces_S1_V5c
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=7
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/traces_S1_V5c_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
S1=/scratch/dhuruva/split1
SR=$S1/refitting_results/fitting/n120/seed20262009
export BCL_EMODEL=delta-split1
export BCL_L5_SIMS=$SR/Sabrina_L5TTPC_L5TTPC_STDP/simulations
export BCL_L5_CIRCUIT=$ROOT/data/dhuruva_delta-split1_circuit_config.json
export BCL_L23_SIMS=$SR/Ebner2019_L23PC_L5TTPC/simulations
export BCL_L23_CIRCUIT=$ROOT/data/dhuruva_split1_l23l5_circuit_config.json
export ANALYTICAL_BASIS_DIR=$S1/basis_l5l5
export L23_BASIS_DIR=$S1/basis_l23l5
OUT=/scratch/dhuruva/bcl_traces_S1_V5c
mkdir -p $OUT
if [ "${STAGE:-traces}" = "plot" ]; then
    python -u glusynapse_v2/bcl_validation/plot_traces_v5.py --out $OUT
else
    python -u glusynapse_v2/bcl_validation/traces_v5.py --check --out $OUT
    python -u glusynapse_v2/bcl_validation/traces_v5.py --out $OUT --workers $((SLURM_CPUS_PER_TASK - 1)) --timeout 1500
fi
