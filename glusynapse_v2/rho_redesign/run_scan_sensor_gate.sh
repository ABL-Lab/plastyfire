#!/bin/bash
#SBATCH --job-name=scan_sensor_gate
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/scan_sensor_gate_%j.out
# CA_DECODE.md section 3: single-pool Ca-sensor gate in place of the C1/C2 c_VDCC gate, no refit, A0g_s3. CPU/numba.
# usage: PILOT=1 sbatch glusynapse_v2/rho_redesign/run_scan_sensor_gate.sh   (Kd 0.5, tau_g 0.5 s; outputs on /scratch)
#        sbatch --time=<from pilot> glusynapse_v2/rho_redesign/run_scan_sensor_gate.sh   (full grid, writes CA_DECODE.md)
# Sizing basis:
#  mem: same L5 load (5 dirs, effcai + vdcc) as scan_vgate_amp 22120612: 19.0 GB MaxRSS -> x 1.25 = 24G. Extras here are
#       per-record (cacr float32 transient, pool decay arrays 4 x T) and the target-row table (< 0.5 GB full grid).
#  time: 22120612 = 3:27 for 2 loads + one replay of 24 configs (12 theta x C1/C2).
#       Full grid = 2 types x 2 Kd x 4 tau_g x 11 theta x C1/C2 = 352 configs (+A0) = 14.7x -> 51 min scaled, plus an
#       extra L2/3 load and pre-pass -> ~55 min -> request 1:30. Too uncertain: run the pilot first.
#       Pilot = 2 types x 1 Kd x 1 tau_g x 11 theta x 2 = 44 configs: 3:27 x 44/24 = 6:20 + L2/3 reload ~1 + pre-pass
#       and window features ~1 -> ~8.5 min x 1.5 = 12.5 -> 0:15. 1 CPU (single-threaded numba, no pool).
#  full run: time = pilot elapsed + 7 x (pilot "L5 replay" + "L2/3 replay" phase seconds, from the log) + 50%.
# MEASURED 22134502 (pilot): TIMEOUT at 15:24 after L5 replay (894 s), 13.91 GB MaxRSS, 91% CPU. Window check answered: win_len alone AUC 0.182, S_kd0.5 mean 0.591 / fix200 0.487 -> integral AUC was a window-length confound. Full grid not run.
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
V2=glusynapse_v2; R=$V2/rho_redesign
FIT=${FIT:-$R/results/v4_A0g_s3.json}
python -m py_compile $R/scan_sensor_gate.py
if [ "${PILOT:-0}" = 1 ]; then
    O=/scratch/dhuruva/sensor_gate_pilot; mkdir -p $O
    python -u $R/scan_sensor_gate.py --fit $FIT --kd 0.5 --taug 500 --save $O/sensor_gate_A0g_s3_pilot \
        --fig $O/fig10_sensor_gate_pilot.png --md $O/CA_DECODE_pilot.md
else
    python -u $R/scan_sensor_gate.py --fit $FIT --save $R/results/sensor_gate_A0g_s3 --fig $R/figs/fig10_sensor_gate.png \
        --md $R/CA_DECODE.md
fi
