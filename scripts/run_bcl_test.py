"""
Standalone bluecellulab simulation matching prefire_simulation_config_test.json.

Applies the same settings as the full edges pipeline:
  - SK_E2 removed (gSK_E2bar = 0)
  - All GluSynapse globals set explicitly after mechanism compilation
  - Per-synapse RNG seeds set to match Neurodamus Random123 seeding

Outputs:
  bcl_test_soma.npy   -- shape (N,2): col0=time[ms], col1=Vsoma[mV]
  bcl_test_rho.npy    -- shape (M,2): col0=time[ms], col1... (one column per synapse)
  bcl_test_rho.txt    -- text summary of initial/final rho per synapse
"""

import json, os
import numpy as np
import bluecellulab
from libsonata import SpikeReader

# Load compiled NMODL mechanisms (CaDynamics_DC0, SK_E2, GluSynapse, etc.)
bluecellulab.neuron.load_mechanisms("/project/ctb-emuller/dhuruva/DEES_cell_packages/")

# ── Config ──────────────────────────────────────────────────────────────────
WORKDIR    = os.path.dirname(os.path.abspath(__file__))
CFG_FILE   = os.path.join(WORKDIR, "prefire_simulation_config_test.json")
SPIKES_H5  = os.path.join(WORKDIR, "prefire_prespikes.h5")
NODE_POP   = "S1nonbarrel_neurons"
PRE_GID    = 180164
POST_GID   = 197248
RECORD_DT  = 0.1   # ms

# ── Load config ──────────────────────────────────────────────────────────────
with open(CFG_FILE) as f:
    cfg = json.load(f)

t_end = cfg["run"]["tstop"]

# Save GluSynapse params BEFORE stripping (needed for manual application after instantiation)
glusyn = dict(cfg.get("conditions", {}).get("mechanisms", {}).get("GluSynapse", {}))

# Strip synapse reports (BCL doesn't support them; we record rho manually)
cfg["reports"] = {k: v for k, v in cfg.get("reports", {}).items()
                  if v.get("type") != "synapse"}

# BCL calls set_global_condition_parameters() before mechanisms are compiled,
# so GluSynapse HOC globals don't exist yet → LookupError for init_depleted/minis_single_vesicle.
# Strip the whole GluSynapse block; we apply all params manually after instantiate_gids().
cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)

import tempfile
tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=WORKDIR, delete=False)
json.dump(cfg, tmp)
tmp.close()

# ── Build circuit ─────────────────────────────────────────────────────────────
BASE_SEED    = cfg["run"].get("random_seed", 0)
SYNAPSE_SEED = cfg["run"].get("synapse_seed", 0)
print(f"Loading circuit (tstop={t_end} ms, base_seed={BASE_SEED}, synapse_seed={SYNAPSE_SEED})...")
sim = bluecellulab.CircuitSimulation(tmp.name, base_seed=BASE_SEED)

# Set synapse_seed explicitly (BCL RNGSettings may not read it from SONATA run block)
from bluecellulab.rngsettings import RNGSettings
rng = RNGSettings.get_instance()
rng.synapse_seed = SYNAPSE_SEED
print(f"RNG: mode={rng.mode}  base_seed={rng.base_seed}  synapse_seed={rng.synapse_seed}")
pre_spikes = SpikeReader(SPIKES_H5)[NODE_POP].get_dict()["timestamps"]
# keep only spikes within simulation window
pre_spikes = pre_spikes[pre_spikes <= t_end]
print(f"Pre-spikes in window: {pre_spikes}")

sim.instantiate_gids(
    [(NODE_POP, POST_GID)],
    add_synapses=True,
    add_minis=False,
    add_pulse_stimuli=True,
    intersect_pre_gids=[(NODE_POP, PRE_GID)],
    pre_spike_trains={(NODE_POP, PRE_GID): pre_spikes},
)
cell = sim.cells[(NODE_POP, POST_GID)]

# ── SK_E2 removal ─────────────────────────────────────────────────────────────
print("Removing SK_E2 (setting gSK_E2bar = 0)...")
for sec in cell.somatic + cell.axonal:
    for seg in sec:
        if hasattr(seg, "gSK_E2bar_SK_E2"):
            seg.gSK_E2bar_SK_E2 = 0.0

# ── Set ALL GluSynapse globals explicitly (after mechanisms are compiled) ──────
# Matches the approach in abl-invivo-plastic-neuron-sscx/ion_cluster_experiments/
# run_spike_pairing_ion.py:90-98. Values from config / fitted chindemi params.
h = bluecellulab.neuron.h
h.cao_CR_GluSynapse             = glusyn.get("cao_CR", 2.0)
h.tau_effca_GB_GluSynapse       = glusyn.get("tau_effca_GB", 278.3177658387)
h.gamma_d_GB_GluSynapse         = glusyn.get("gamma_d_GB", 101.5387594661)
h.gamma_p_GB_GluSynapse         = glusyn.get("gamma_p_GB", 216.1841700668)
h.init_depleted_GluSynapse      = float(glusyn.get("init_depleted", True))
h.minis_single_vesicle_GluSynapse = float(glusyn.get("minis_single_vesicle", False))
print(f"  cao_CR_GluSynapse             = {h.cao_CR_GluSynapse}")
print(f"  tau_effca_GB_GluSynapse       = {h.tau_effca_GB_GluSynapse}")
print(f"  gamma_d_GB_GluSynapse         = {h.gamma_d_GB_GluSynapse}")
print(f"  gamma_p_GB_GluSynapse         = {h.gamma_p_GB_GluSynapse}")
print(f"  init_depleted_GluSynapse      = {h.init_depleted_GluSynapse}")
print(f"  minis_single_vesicle_GluSynapse = {h.minis_single_vesicle_GluSynapse}")

# ── Per-synapse RNG seeding to match Neurodamus Random123 ─────────────────────
# Mirrors abl-invivo-plastic-neuron-sscx/ion_cluster_experiments/
# run_spike_pairing_ion.py:209-214.
# ND seeds: (tgid=POST_GID, 100000+synapse_local_idx, synapse_seed+200)
print("\nSetting per-synapse RNG seeds...")
for syn_id, synapse in cell.synapses.items():
    local_idx = int(syn_id[1])
    s1 = POST_GID
    s2 = 100000 + local_idx
    s3 = 200 + SYNAPSE_SEED
    synapse.randseed1 = s1
    synapse.randseed2 = s2
    synapse.randseed3 = s3
    synapse.hsynapse.setRNG(s1, s2, s3)
    print(f"  syn {syn_id}: setRNG({s1}, {s2}, {s3})")

# ── Recording ─────────────────────────────────────────────────────────────────
# Soma voltage
soma_v  = bluecellulab.neuron.h.Vector()
soma_t  = bluecellulab.neuron.h.Vector()
soma_v.record(cell.soma(0.5)._ref_v)
soma_t.record(bluecellulab.neuron.h._ref_t)

# Rho per synapse
rho_vecs = {}
for syn_id, synapse in cell.synapses.items():
    v = bluecellulab.neuron.h.Vector()
    v.record(synapse.hsynapse._ref_rho_GB)
    rho_vecs[syn_id] = v

# Print initial state
print("\nInitial synapse state:")
for syn_id, syn in cell.synapses.items():
    h = syn.hsynapse
    print(f"  syn {syn_id}: rho0={h.rho0_GB:.3f}  theta_d={h.theta_d_GB:.4f}  theta_p={h.theta_p_GB:.4f}")

# ── Run ───────────────────────────────────────────────────────────────────────
print(f"\nRunning BCL simulation for {t_end} ms (CVODE)...")
sim.run(t_end, cvode=True)
print("Done.")

# ── Save ──────────────────────────────────────────────────────────────────────
t_arr   = np.array(soma_t)
v_arr   = np.array(soma_v)
soma_out = np.column_stack([t_arr, v_arr])
np.save(os.path.join(WORKDIR, "bcl_test_soma.npy"), soma_out)
print(f"Saved bcl_test_soma.npy  shape={soma_out.shape}")

# Rho traces (interpolate to common time axis)
t_common = np.arange(0, t_end + RECORD_DT, RECORD_DT)
rho_out  = np.zeros((len(t_common), 1 + len(rho_vecs)))
rho_out[:, 0] = t_common
for i, (syn_id, v) in enumerate(rho_vecs.items()):
    rho_out[:, i+1] = np.interp(t_common, t_arr, np.array(v))
np.save(os.path.join(WORKDIR, "bcl_test_rho.npy"), rho_out)
print(f"Saved bcl_test_rho.npy   shape={rho_out.shape}")

# Text summary
with open(os.path.join(WORKDIR, "bcl_test_rho.txt"), "w") as fh:
    fh.write(f"{'syn_id':<20}  {'initial_rho':>12}  {'final_rho':>10}\n")
    for i, (syn_id, v) in enumerate(rho_vecs.items()):
        rv = np.array(v)
        fh.write(f"{str(syn_id):<20}  {rv[0]:>12.6f}  {rv[-1]:>10.6f}\n")
print("Saved bcl_test_rho.txt")

# Cleanup temp config
os.unlink(tmp.name)
print("\nDone. Compare bcl_test_soma.npy with ND soma report to see timing delta.")
