#!/bin/bash
# NOT SUBMITTED (est. ~24 CPU h; part A + B ~72 > 30, orchestrator logs first). t_drive 3 full re-prefire, part B: Sjostrom 2001 (10 protocols) +
# Markram (7 protocols), 24 subset24 pairs, hash delta-prefire-vseg. 3 array tasks x 8 pairs x 17 protocols = 136 sims.
# MEASURED 22073523: MaxRSS 103-110 GiB (at the 110G limit!), 18-25 min, 74% CPU eff on 16 -> next run 140G, 0:45.
# Sizing (no v_seg measurement for these protocols; prior per-sim times x1.32 v_seg factor measured in pilot 22073044):
#   Markram 105 s (21903173) -> 136 s; sj01 0.1 Hz +-10 ~375 s (logs) -> 490 s; sj01 10-50 Hz ~110 s -> 143 s.
#   per pair 955 + 2120 = 3075 s; 8 pairs / 16 CPUs = 26 min + tail = ~30 min. Time 30 + 50% = 0:45.
#   Memory 5.5 GB/worker (21903173: 240 GB / 44) x 16 = 88 GB + 25% = 110G.
#SBATCH --job-name=vseg_sj01mk
#SBATCH --account=rrg-emuller
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=140G
#SBATCH --time=00:45:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vseg_sj01mk_%A_%a.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(python -c "
p=open('glusynapse_v2/subset24_pairs.txt').read().strip().split(',')
print(','.join(p[${SLURM_ARRAY_TASK_ID:-0}::3]))")
PROTOS=10Hz_-50ms,10Hz_-30ms,10Hz_-10ms,10Hz_5ms,10Hz_10ms,10Hz_30ms,10Hz_50ms,sjostrom_0.1hz_dt+10ms,sjostrom_0.1hz_dt-10ms,sjostrom_10hz_dt+10ms,sjostrom_10hz_dt-10ms,sjostrom_20hz_dt+10ms,sjostrom_20hz_dt-10ms,sjostrom_40hz_dt+10ms,sjostrom_40hz_dt-10ms,sjostrom_50hz_dt+10ms,sjostrom_50hz_dt-10ms
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 16 --protocols $PROTOS
