"""Config was ruled out (A/B/C all found 1.60 nA). Remaining variable: the
c_pre simulation that _find_cpre_cpost runs in the SAME process before the
c_post threshold search.

Bisects three states, each in a fresh process (run one stage per invocation):
  stage0: search alone                         -> known to PASS
  stage1: search after instantiate_gids only   (no run)
  stage2: search after instantiate + run(1500) (full c_pre replication)
  stage3: search after the PRODUCTION sim is instantiated
          (add_pulse_stimuli=True + full prespike replay) -- what the real
          _run_prefire_edges_process has standing before it calls _find_cpre_cpost

Usage: python diag_cpost_state.py <stage> [gid] [pair]
"""
import json
import os
import sys
import tempfile

import bluecellulab
from libsonata import SpikeReader
from plastyfire.simulator import spike_threshold_finder, _set_global_params

STAGE = int(sys.argv[1]) if len(sys.argv) > 1 else 0
GID = int(sys.argv[2]) if len(sys.argv) > 2 else 205995
PAIR = sys.argv[3] if len(sys.argv) > 3 else "13255-205995"
PRE_GID = int(PAIR.split("-")[0])

BASE = ("/lustre06/project/6077694/dhuruva/plastyfire/refitting_results/"
        "fitting/n100/seed19091997/L23PC_L5TTPC_STDP/simulations")
PROD = os.path.join(BASE, PAIR, "10Hz_5ms", "prefire_simulation_config.json")
NODE_POP = "S1nonbarrel_neurons"

FIT_PARAMS = {"gamma_d_GB_GluSynapse": 77.7558, "gamma_p_GB_GluSynapse": 299.9121,
              "tau_effca_GB_GluSynapse": 278.3177658387}

# same stripping _find_cpre_cpost does before handing the config to the search
cfg = json.load(open(PROD))
cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
cfg.pop("reports", None)
cfg.setdefault("output", {})["output_dir"] = tempfile.mkdtemp()
tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
json.dump(cfg, tmp)
tmp.close()
sim_config_eff = tmp.name
seed = cfg["run"]["random_seed"]

print(f"stage {STAGE}, gid {GID}, pair {PAIR}", flush=True)

if STAGE in (1, 2):
    print("  building c_pre sim (instantiate_gids)…", flush=True)
    sim_pre = bluecellulab.CircuitSimulation(sim_config_eff, base_seed=seed)
    sim_pre.instantiate_gids(
        [(NODE_POP, GID)],
        add_synapses=True, add_minis=False, add_pulse_stimuli=False,
        intersect_pre_gids=[(NODE_POP, PRE_GID)],
        pre_spike_trains={(NODE_POP, PRE_GID): [1000.0]},
    )
    cell_pre = sim_pre.cells[(NODE_POP, GID)]
    for sec in cell_pre.somatic + cell_pre.axonal:
        sec.uninsert("SK_E2")
    _set_global_params(FIT_PARAMS)
    print(f"  instantiated, {len(cell_pre.synapses)} synapses", flush=True)

if STAGE == 2:
    print("  running c_pre sim to t=1500…", flush=True)
    sim_pre.run(1500.0, cvode=True)
    print(f"  done, h.t = {bluecellulab.neuron.h.t}", flush=True)

if STAGE >= 3:
    # faithful replication of _run_prefire_edges_process: the FULL production sim
    # (induction pulse train + complete prespike replay) is instantiated before
    # _find_cpre_cpost is ever called.
    print("  building PRODUCTION sim (add_pulse_stimuli=True, full replay)…", flush=True)
    pre_spikes_path = os.path.join(BASE, PAIR, "10Hz_5ms", "prefire_prespikes.h5")
    pre_spikes = SpikeReader(pre_spikes_path)[NODE_POP].get_dict()["timestamps"]
    sim_prod = bluecellulab.CircuitSimulation(sim_config_eff, base_seed=seed)
    sim_prod.instantiate_gids(
        [(NODE_POP, GID)],
        add_synapses=True, add_minis=False, add_pulse_stimuli=True,
        intersect_pre_gids=[(NODE_POP, PRE_GID)],
        pre_spike_trains={(NODE_POP, PRE_GID): pre_spikes},
    )
    cell_prod = sim_prod.cells[(NODE_POP, GID)]
    for sec in cell_prod.somatic + cell_prod.axonal:
        sec.uninsert("SK_E2")
    _set_global_params(FIT_PARAMS)
    print(f"  production cell instantiated, {len(cell_prod.synapses)} synapses, "
          f"{len(pre_spikes)} prespikes", flush=True)

print("  now running the c_post single-AP search…", flush=True)
for w in [1.5, 3, 5]:
    res = spike_threshold_finder(sim_config_eff, GID, 1, 0.1, w, 1000.,
                                 0.05, 5., 100, NODE_POP, True)
    if res is not None:
        print(f"  -> width {w} ms: FOUND amp={res['amp']:.4f} nA, "
              f"{len(res['t_spikes'])} spike(s)", flush=True)
        break
    print(f"  -> width {w} ms: no single AP up to 5.0 nA", flush=True)
else:
    print("  -> FAILED at every width (reproduced the bug)", flush=True)
