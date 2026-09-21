# Edges-Based Fitting Pipeline

**Created:** 2026-06-04  
**Repo:** `/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire`

---

## Overview

Three new files that replicate `modelfitter.py` / `evaluator.py` but drive induction simulations via `dhuruva_modified_edges.h5` (bluecellulab prefire pipeline) instead of the recipe-based EPG pipeline.

---

## Created Files

### `plastyfire/pairrunner_edges_fit.py`

Thin CLI wrapper analogous to `pairrunner.py` for the edges pipeline.

- Accepts all 10 fitted params as individual `--name=value` args
- Calls `simulator_edges.runconnectedpair_prefire_from_edges()`
- Saves output to `simulation_edges_{param_hash}.pkl` in the workdir
- **Key role:** subprocess target called by the Pool worker — as a fresh OS process (not a daemon) it can freely spawn the internal `multiprocessing.Process` used by the simulator, avoiding nested-daemon errors

```
python pairrunner_edges_fit.py \
    --gamma_d_GB_GluSynapse=120.0 --gamma_p_GB_GluSynapse=210.0 \
    --a00=2.0 --a01=1.5 --a10=2.0 --a11=2.0 --a20=3.0 --a21=2.0 --a30=4.0 --a31=2.0 \
    --fastforward=280000 \
    --edges-h5 /path/to/dhuruva_modified_edges.h5 \
    --param_hash abc123def456
```

---

### `plastyfire/evaluator_edges.py`

Contains three public symbols:

#### `EDGES_FIT_PARAM_NAMES`
The 10 parameters optimised (tau_effca is fixed at `FITTED_TAU = 278.3177658387`):
```
gamma_d_GB_GluSynapse, gamma_p_GB_GluSynapse,
a00, a01, a10, a11, a20, a21, a30, a31
```

#### `run_simulation_worker_edges(args)`
Module-level Pool worker function (must be top-level for pickle).  
`args` tuple: `(param_values, sim_dict, param_hash, edges_h5, cpre_cpost_cache)`  
Dispatches `pairrunner_edges_fit.py` via `subprocess.run()`, then calls `compute_epsp_ratio_edges_batch`.

#### `compute_epsp_ratio_edges_batch(param_values, sim_dict, workdir)`
- Reads `simulation_edges_{hash}.pkl` (has flat `initial_rho` / `final_rho` float lists)
- Thresholds rho values at 0.5 → binary state
- Looks up pre-existing `ephys_data_{pre}_{post}_{rho_str}.pkl` files 4 levels above workdir
- Returns `(protocol_id, epsp_ratio, [])` or `(None, None, [])` on failure
- **Does NOT auto-generate ephys files** — see prerequisite note below

#### `EvaluatorEdges(Evaluator)`
Inherits all setup logic from `Evaluator` (`__init__`, objectives, `all_sims`, scoring).  
Adds constructor args: `edges_h5`, `cpre_cpost_cache`, `basis_dir`.  
Overrides: `evaluate_with_multiprocessing`.

---

### `plastyfire/modelfitter_edges.py`

Mirror of `modelfitter.py` using `EvaluatorEdges`.

**Additional CLI args (vs `modelfitter.py`):**

| Arg | Default | Description |
|-----|---------|-------------|
| `--edges-h5` | `EDGES_H5_DEFAULT` | Path to `dhuruva_modified_edges.h5` |
| `--cpre-cpost-cache` | `None` | Precomputed c_pre/c_post cache pkl |
| `--basis-dir` | **required** | Directory with `basis_{pre}_{post}.csv` files from `run_basis_pair_edges.py` |

**Output files:** `checkpoint_edges.pkl`, `bestsol_edges.pkl` (instead of `checkpoint.pkl`, `bestsol.pkl`)

**Example run:**
```bash
python plastyfire/modelfitter_edges.py \
    --edges-h5 /project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5 \
    --basis-dir basis_results_edges_mini \
    --cpre-cpost-cache /path/to/cpre_cpost_cache.pkl \
    --sample_size 100 --pop_size 128 --gen 30 \
    --use-multiprocessing
```

---

## Design Decisions & Known Issues

### 1. Nested Multiprocessing (solved)
`simulator_edges.runconnectedpair_prefire_from_edges()` internally spawns a `multiprocessing.Process`. Pool workers are daemon processes and cannot fork. **Fix:** Pool workers call `pairrunner_edges_fit.py` via `subprocess.run()` — the OS-level subprocess is not a daemon and can fork freely.

### 2. EPSP Ratio — Basis Technique (no ephys files needed)

Instead of running bluecellulab test-pulse simulations, the EPSP ratio is computed analytically from pre-generated basis CSVs.

**`run_basis_pair_edges.py`** (already run → `basis_results_edges_mini/`) measures EPSP for:
- all-depressed config: `"0,0,...,0"` → **e₀**
- each singleton: one rho_i=1, rest=0 → **eᵢ_singleton**
- all-potentiated: `"1,1,...,1"` (validation)

During fitting, the EPSP for any rho vector is computed analytically:
```
EPSP(ρ) = e₀ + Σᵢ [ ρᵢ × (eᵢ_singleton - e₀) ]
```
This uses the **continuous rho values** directly from the induction simulation (no 0.5 threshold), since gmax scales linearly with rho and EPSPs sum linearly across synapses.

**Result:** No `ephys_data_*.pkl` files needed at all. EPSP ratio is pure numpy arithmetic — microseconds per evaluation.

> **Note:** Linear superposition slightly overestimates EPSP for high-rho configs (~12% at all-ones) due to driving-force saturation, but this is the same approximation used by `plot_refitting_stdp_basis.py` and is systematic across all candidates.

### 3. Bug in original `evaluator.py` (not fixed, documented)
`run_simulation_worker` uses module-level `FIT_PARAM_NAMES` (11 items, starts with `tau_effca`) to map `param_values`, but `modelfitter.py`'s `FIT_PARAMS` only has 10 params (tau_effca commented out). This causes all param assignments to be offset by one. The edges evaluator avoids this bug by using the correct `EDGES_FIT_PARAM_NAMES` (10 items, no tau_effca) and always injecting `tau_effca=FITTED_TAU` explicitly.

---

## Pipeline Comparison

| Step | Original (`modelfitter.py`) | Edges (`modelfitter_edges.py`) |
|------|----------------------------|-------------------------------|
| Simulation script | `pairrunner.py` | `pairrunner_edges_fit.py` |
| Simulation engine | `simulator.runconnectedpair_induction()` | `simulator_edges.runconnectedpair_prefire_from_edges()` |
| Config file | `simulation_config.json` | `prefire_simulation_config.json` |
| Synapse params source | EPG recipe calibration | `dhuruva_modified_edges.h5` (auto-loaded by bluecellulab) |
| Threshold injection | EPG recipe | From `edges.h5` directly or computed via a-params |
| Result pkl format | `rho_GB` time-series matrix | `initial_rho` / `final_rho` flat float lists |
| Result pkl name | `simulation_{hash}.pkl` | `simulation_edges_{hash}.pkl` |
| EPSP ratio source | `ephys_data_*.pkl` + bluecellulab | Basis CSV + linear superposition (no simulation) |
| Checkpoint file | `checkpoint.pkl` | `checkpoint_edges.pkl` |
| Best solution file | `bestsol.pkl` | `bestsol_edges.pkl` |
