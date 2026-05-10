"""
Diagnostic: print synapse ordering from bluecellulab cell vs basis CSV.
Run in the bluecellulab cluster environment from the plastyfire root.

Usage:
    python diagnose_syn_order.py
"""
import os, sys, json
import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, "/lustre06/project/6077694/dhuruva/plastyfire")
import bluecellulab

MECHANISMS_PATH = "/project/ctb-emuller/dhuruva/DEES_cell_packages/"
bluecellulab.neuron.load_mechanisms(MECHANISMS_PATH)

PAIR         = "180164-197248"
DT           = "10Hz_5ms"
PRE_GID      = 180164
POST_GID     = 197248
NODE_POP     = "S1nonbarrel_neurons"
EDGE_POP     = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
EDGES_H5     = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
WORKDIR      = (f"/lustre06/project/6077694/dhuruva/plastyfire/refitting_results/"
                f"fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations/{PAIR}/{DT}")
BASIS_CSV    = (f"/lustre06/project/6077694/dhuruva/plastyfire/"
                f"basis_results_edges_mini/basis_180164_197248.csv")

# ── 1. Load BCL cell ──────────────────────────────────────────────────────────
sim_config = os.path.join(WORKDIR, "prefire_simulation_config.json")
with open(sim_config) as f:
    cfg = json.load(f)
cfg.pop("reports", None)
cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
import tempfile
tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=WORKDIR, delete=False)
json.dump(cfg, tmp); tmp.close()

sim = bluecellulab.CircuitSimulation(tmp.name)
sim.instantiate_gids(
    [(NODE_POP, POST_GID)],
    add_synapses=True, add_minis=False, add_pulse_stimuli=False,
    intersect_pre_gids=[(NODE_POP, PRE_GID)],
)
os.unlink(tmp.name)
cell = sim.cells[(NODE_POP, POST_GID)]

# ── 2. Print BCL synapse order: position → (local_id, gmax_d, gmax_p, rho0) ──
print("\n=== BCL cell.synapses iteration order ===")
print(f"{'pos':>4}  {'syn_id':>8}  {'gmax_d':>8}  {'gmax_p':>8}  {'rho0':>6}")
bcl_order = []
for pos, (syn_id, synapse) in enumerate(cell.synapses.items()):
    h = synapse.hsynapse
    gmax_d = float(h.gmax_d_AMPA)
    gmax_p = float(h.gmax_p_AMPA)
    rho0   = float(h.rho0_GB)
    print(f"{pos:>4}  {str(syn_id):>8}  {gmax_d:>8.4f}  {gmax_p:>8.4f}  {rho0:>6.1f}")
    bcl_order.append({"pos": pos, "syn_id": syn_id, "gmax_d": gmax_d, "gmax_p": gmax_p, "rho0": rho0})

# ── 3. Load global_ids from rho.h5 ───────────────────────────────────────────
rho_h5 = os.path.join(WORKDIR, "bluecellulab_results", "rho.h5")
if os.path.exists(rho_h5):
    with h5py.File(rho_h5) as f:
        eids = f[f"report/{NODE_POP}/mapping/element_ids"][()]
        data = f[f"report/{NODE_POP}/data"][()]
    print(f"\n=== rho.h5 element_ids (global circuit IDs) ===")
    print(f"{'pos':>4}  {'element_id':>12}  {'initial_rho':>11}  {'final_rho':>10}")
    for pos, (eid, r0, r1) in enumerate(zip(eids, data[0], data[-1])):
        print(f"{pos:>4}  {eid:>12}  {r0:>11.4f}  {r1:>10.4f}")

    print(f"\n=== BCL cell.synapses pos vs rho.h5 element_ids: do they correspond? ===")
    with h5py.File(EDGES_H5) as f:
        pop0 = f"edges/{EDGE_POP}/0"
        gmax_d_arr = f[f"{pop0}/gmax_d_AMPA"][eids.tolist()]
        gmax_p_arr = f[f"{pop0}/gmax_p_AMPA"][eids.tolist()]
        rho0_arr   = f[f"{pop0}/rho0_GB"][eids.tolist()]
    print(f"{'pos':>4}  {'element_id':>12}  {'edges gmax_p':>13}  {'bcl gmax_p':>10}  {'match':>6}")
    for pos, (eid, gp_edges, bcl) in enumerate(zip(eids, gmax_p_arr, bcl_order)):
        match = abs(gp_edges - bcl['gmax_p']) < 1e-3
        print(f"{pos:>4}  {eid:>12}  {gp_edges:>13.4f}  {bcl['gmax_p']:>10.4f}  {'OK' if match else 'MISMATCH':>6}")

# ── 4. Basis singleton EPSPs ──────────────────────────────────────────────────
if os.path.exists(BASIS_CSV):
    df = pd.read_csv(BASIS_CSV)
    df['cl'] = df['config'].apply(lambda x: [int(i) for i in x.split(',')])
    print(f"\n=== Basis singletons: which position has which EPSP ===")
    print(f"{'pos':>4}  {'singleton_mean':>15}")
    for _, row in df[df['cl'].apply(lambda x: sum(x)==1)].sort_values('config').iterrows():
        pos = row['cl'].index(1)
        print(f"{pos:>4}  {row['mean']:>15.4f}")

# ── 5. Repeat synapse loading with simulation_config.json (used by basis) ────
SIM_CONFIG_BASIS = os.path.join(WORKDIR, "simulation_config.json")
if os.path.exists(SIM_CONFIG_BASIS):
    print(f"\n=== BCL cell.synapses order using simulation_config.json (basis config) ===")
    with open(SIM_CONFIG_BASIS) as f:
        cfg2 = json.load(f)
    cfg2.pop("reports", None)
    cfg2.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
    tmp2 = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=WORKDIR, delete=False)
    json.dump(cfg2, tmp2); tmp2.close()

    sim2 = bluecellulab.CircuitSimulation(tmp2.name)
    sim2.instantiate_gids(
        [(NODE_POP, POST_GID)],
        add_synapses=True, add_minis=False, add_pulse_stimuli=False,
        intersect_pre_gids=[(NODE_POP, PRE_GID)],
    )
    os.unlink(tmp2.name)
    cell2 = sim2.cells[(NODE_POP, POST_GID)]

    print(f"{'pos':>4}  {'syn_id':>8}  {'gmax_d':>8}  {'gmax_p':>8}  {'rho0':>6}  {'matches_prefire':>15}")
    for pos, (syn_id, synapse) in enumerate(cell2.synapses.items()):
        h2 = synapse.hsynapse
        gmax_p_basis = float(h2.gmax_p_AMPA)
        match = abs(gmax_p_basis - bcl_order[pos]['gmax_p']) < 1e-3 if pos < len(bcl_order) else False
        print(f"{pos:>4}  {str(syn_id):>8}  {float(h2.gmax_d_AMPA):>8.4f}  {gmax_p_basis:>8.4f}  "
              f"{float(h2.rho0_GB):>6.1f}  {'OK' if match else 'MISMATCH':>15}")
