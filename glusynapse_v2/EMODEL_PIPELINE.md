# Re-running the plasticity pipeline on a new L5 TTPC emodel

One command (details in the header of `run_new_emodel_chain.sh`):

    EMODEL_NAME=<name> EMODEL_DIR=<dir> [CFG_L5L5/CFG_L23L5/CFG_L23L23] [PATHWAYS=l5l5,l23l5,l23l23] [PILOT=1] [SKIPFIT=1] [DRYRUN=1] [DECISIONS_LOGGED=1] bash glusynapse_v2/run_new_emodel_chain.sh

Prerequisite: the new mechanisms (new .mod, new SUFFIX) must be compiled into `/project/rrg-emuller/dhuruva/DEES_cell_packages/x86_64/libnrnmech`
(path hard-coded in `plastyfire/simwriter.py:28` and `simulator_edges.py:37`). Raw output goes to `/scratch/dhuruva/<name>/`, npz to
`extracted/<dataset>_delta-<name>-prefire-{vca|vseg-rs}/`. og-delta outputs are never touched (reserved names and the og dir are refused).

## Where the emodel is referenced (all of it)
1. `biophysical_neuron_models_dir` in the circuit config (`data/dhuruva_delta_circuit_config.json`, `_l23l5_`, `_l23l23_`; they differ only in the edges file).
   The driver writes `data/dhuruva_<name>[_l23l5|_l23l23]_circuit_config.json` by replacing only this key. Edges files (`dhuruva_modified_edges*.h5`) are emodel-independent.
2. The compiled mechanisms (point 1 of the prerequisite).
3. The simwriter yaml `circuit.config` (also the calibration config); the driver rewrites it in derived yamls.
Nothing else: the param hash (`delta-<name>-prefire-vseg`) is only a label; fit scripts take directories.

## Steps (sizes are measured og-delta equivalents; CPU·h = est. for the full run)
Emodel-dependent = y means the output changes with the emodel. The L5 set uses the 24-pair Sjostrom2007 subset; the EPSP basis is per all pairs (120 for L2/3->L5).

| Pathway | Step | Script | Inputs | Emodel dep. | Sizing (cpu / mem / time) | CPU·h |
|---|---|---|---|---|---|---|
| L5->L5 Sabrina/Markram | sim_l5_mk | plastyfire/simwriter.py | yaml, pairs, circuit cfg | y (threshold calibration, single_cells pkl) | 8 / 2G / 0:15 | 0.3 |
| L5->L5 Ebner (sj01) | sim_l5_sj01 | simwriter.py | yaml | y | 8 / 2G / 0:30 | 0.8 |
| L5->L5 Sjostrom03 (+r50) | sim_l5_sj03, sj03r50 | simwriter.py | yaml | y | 8 / 2G / 0:30 | 1.6 |
| L5->L5 Sjostrom07 | sim_l5_sj07 | simwriter.py | yaml | y | 8 / 2G / 0:15 (22086838) | 0.7 |
| L5->L5 Sjostrom04 (optional) | sim_l5_sj04 | simwriter.py | yaml | y | 8 / 2G / 0:15 (22127910) | 0.3 |
| L5->L5 | cache_l5 | precompute_cpre_cpost.py | mk workdirs, edges | y (own single-AP search) | 8 / 13G / 0:15 | 0.3 |
| L5->L5 | basis_l5 | run_basis_pair_edges.py (array per pair) | mk workdirs | y | 12 / 15G / 0:30 per pair | 29 |
| L5->L5 | pf_l5_A (mk+sj01) | run_de_fit2_pool.py --cooker (array 3) | sims, cache | y | 16 / 140G / 0:45 | 19 |
| L5->L5 | pf_l5_B (sj03+r50) | same | same | y | 32 / 205G / 0:45 | 37.6 |
| L5->L5 | pf_l5_sj07, pf_l5_sj04 | same | same | y | 16 / 98G / 0:45; 12 / 70G / 0:30 | 5.8; 2.3 |
| L5->L5 | ext_l5 | extract_v2.py | prefire traces, cache | y | 4 / 145G / 0:30 (22126528, 22074074) | 1.5 |
| L2/3->L5 | sim_l23l5 (+nulls) | simwriter.py | yaml, pairs from index csv | y | 16 / 25G / 1:00 (22043887); nulls 1 / 1G | 8 |
| L2/3->L5 | cache_l23l5 | precompute_cpre_cpost.py | workdirs | y | 16 / 25G / 0:15 (22085683) | 1.4 |
| L2/3->L5 | basis_l23l5 | run_basis_pair_edges.py | workdirs | y | 12 / 15G / 0:30 per pair | 170 |
| L2/3->L5 | pf_l23l5 (+nulls) | run_de_fit2_pool.py | sims, cache | y | 32 / 255G / 1:00 (22089155); nulls 16 / 150G / 1:00 | 20; 10.5 |
| L2/3->L5 | ext_l23l5 | extract_v2.py | traces | y | 4 / 99G / 0:15 | 0.6 |
| L2/3->L2/3 | sim_l23l23 | simwriter.py | yaml, pairs | y | 16 / 25G / 0:30 (22128319) | 5.1 |
| L2/3->L2/3 | cache_l23l23 | precompute_cpre_cpost.py | workdirs | y | 16 / 13G / 0:15 (22128321) | 0.7 |
| L2/3->L2/3 | basis_l23l23 | run_basis_pair_edges.py | workdirs | y | 12 / 7G / 0:15 per pair (22128322) | 24 |
| L2/3->L2/3 | pf_l23l23 | run_de_fit2_pool.py | sims, cache | y | 32 / 310G / 1:45 | 39.4 |
| L2/3->L2/3 | ext_l23l23 | extract_v2.py | traces | y | 4 / 124G / 0:45 | 1.4 |
| Joint | fit_s5 (seeded), fit_s6 (unseeded) | rho_redesign/run_fit_v4.sh (fit_v4n.py) | npz of all 3 pathways, basis dirs | no (consumes the above) | GPU 3g.40gb, 1 CPU / 79G / 1:45 (22134473) each | ~1 each |

Total, full run (3 pathways, mk+sj01+sj03+sj07, fit included): about 379 CPU·h (driver estimate), of which the two EPSP bases are 223, the three prefire stages ~135.
Pilot (PILOT=1, 2 pairs x 2 protocols per pathway): about 7 CPU·h for l5l5+l23l5. Steps that do not depend on the emodel: only the fit itself and the yaml/protocol definitions.

## Findings and caveats
- The og-delta Markram (10 Hz) workdirs were calibrated on `data/dhuruva_modified_ion_channels_circuit_config.json` (an older emodel), not on delta. A bit-identical
  reproduction of the Markram npz therefore needs `L5_MK_CAL_CFG=$PWD/data/dhuruva_modified_ion_channels_circuit_config.json`; by default the driver calibrates on the NEW emodel (what a new emodel needs).
- The fit step uses the then-best settings of C1Ajn_s5 (FIT_SET/FIT_FREE/FIT_DROPT/FIT_SEEDFITS overrideable) plus an unseeded basin check.
- Full runs over 30 CPU·h are refused unless `DECISIONS_LOGGED=1` (log them in DECISIONS.md first).
- Plumbing check: with name `ogtest*` and the og-delta dir the driver adds chained compare steps (`compare_emodel_npz.py`, 1 CPU, 2G, 0:15) that print per-record rho/ratios against the existing og-delta extraction.
