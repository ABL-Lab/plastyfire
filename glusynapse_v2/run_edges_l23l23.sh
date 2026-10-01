#!/bin/bash
# MEASURED 22128320: 3:12, MaxRSS 29.5 GB; CHECK OK (296 syn / 117 pairs changed, all L23L23); 120 pairs, 489 syn, rho0=1 0.544 -> next 0:15.
# NEW edges file data/dhuruva_modified_edges_l23l23.h5: the same per-connection conductance median split (k_u 0.2,
# k_gsyn 2, pot 0.5) as dhuruva_modified_edges_l23l5.h5 (T29), on L5-L5 (Sabrina n120) + L2/3->L5 (Ebner L23) +
# L2/3->L2/3 (Zilberter) pairs. dhuruva_modified_edges*.h5 are not touched. Config data/dhuruva_delta_l23l23_circuit_config.json.
# Then a check: vs the l23l5 file, only synapses of L2/3->L2/3 pairs may differ; rho0 / conductance stats of those pairs.
# Sizing from create_dhuruva_edges 21696388 / 21703522 (pair mode, 2 CPUs, 39% eff): 4:37 / 2:24, MaxRSS 24.3 / 30.5 GiB
# (5 full patched arrays, ~9.8 GB file) -> 1 CPU, 39G (30.5 + 25%); 0:30 (create 7 min + check, unmeasured, reads 2 x 9.8 GB).
#SBATCH --job-name=edges_l23l23
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=39G
#SBATCH --time=00:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/edges_l23l23_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=$ROOT/refitting_results/fitting/n120/seed20262009
python -u create_dhuruva_edges.py --output data/dhuruva_modified_edges_l23l23.h5 --no-circuit-config --no-patch-configs \
    --extra-sims-dir $S/Ebner2019_L23PC_L5TTPC/simulations --extra-sims-dir $S/Zilberter2009_L23PC_L23PC/simulations \
    | { grep -v "^  pair \|^    [a-zA-Z0-9_]*: [0-9]*/" || true; }
python -u - <<'EOF'
import os, h5py, numpy as np
R = '/project/rrg-emuller/dhuruva/plastyfire'; S = f'{R}/refitting_results/fitting/n120/seed20262009'
pop = 'edges/S1nonbarrel_neurons__S1nonbarrel_neurons__chemical'; P = pop + '/0/'
F = ['rho0_GB', 'Use_d_TM', 'Use_p_TM', 'gmax_d_AMPA', 'gmax_p_AMPA']
a = h5py.File(f'{R}/data/dhuruva_modified_edges_l23l5.h5', 'r'); b = h5py.File(f'{R}/data/dhuruva_modified_edges_l23l23.h5', 'r')
n = a[P + 'rho0_GB'].shape[0]; diff = np.zeros(n, bool); C = 50_000_000
for f in F:
    for s in range(0, n, C):
        diff[s:s + C] |= a[P + f][s:s + C] != b[P + f][s:s + C]
idx = np.where(diff)[0]
src, tgt = b[pop + '/source_node_id'], b[pop + '/target_node_id']
ch = set(zip(src[idx].tolist(), tgt[idx].tolist())) if len(idx) else set()
pairs = lambda d: {tuple(map(int, x.split('-'))) for x in os.listdir(f'{S}/{d}/simulations') if '-' in x}
l23 = pairs('Zilberter2009_L23PC_L23PC'); old = pairs('Sabrina_L5TTPC_L5TTPC_STDP') | pairs('Ebner2019_L23PC_L5TTPC')
print(f'changed synapses {len(idx)}, changed pairs {len(ch)}, in L23L23 set {len(ch & l23)}, in L5L5/L23L5 sets {len(ch & old)}')
ok = len(ch) > 0 and not (ch & old) and ch <= l23
print('CHECK', 'OK' if ok else 'FAIL')
if not ok:
    raise SystemExit(1)   # stops the afterok chain (cache, basis, pilot)
ix = b[pop + '/indices/target_to_source']; n2r = ix['node_id_to_ranges'][:]; r2e = ix['range_to_edge_id'][:]
cond_ds = b[P + 'conductance']; rho_ds = b[P + 'rho0_GB']; st_ds = b[P + 'syn_type_id']
rho, cond = [], []
for pre, post in sorted(l23):
    for e0, e1 in r2e[n2r[post][0]:n2r[post][1]]:
        e0, e1 = int(e0), int(e1); m = (src[e0:e1] == pre) & (st_ds[e0:e1] >= 100)
        rho.append(rho_ds[e0:e1][m]); cond.append(cond_ds[e0:e1][m])
rho, cond = np.concatenate(rho), np.concatenate(cond)
print(f'L23L23: {len(l23)} pairs, {len(rho)} synapses, rho0=1 fraction {rho.mean():.3f}, '
      f'mean conductance rho0=1 {cond[rho == 1].mean():.3f} vs rho0=0 {cond[rho == 0].mean():.3f} nS')
EOF
