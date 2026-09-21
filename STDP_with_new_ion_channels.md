# STDP Simulations with Modified Ion Channel Models

## Goal

Re-run the BCL STDP induction pipeline (chindemi params) using a modified circuit config
that swaps the neuron biophysical models from the default `emodels_hoc` to
`modified_emodels_hoc` — without disturbing any existing simulation results or the
default code paths.

## What changed in the circuit config

**Old:** `/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_circuit_config.json`
```
"biophysical_neuron_models_dir": "$BASE_DIR/emodels_hoc"
```

**New:** `/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_ion_channels_circuit_config.json`
```
"biophysical_neuron_models_dir": "/project/ctb-emuller/dhuruva/DEES_cell_packages/modified_emodels_hoc"
```

That is the **only** difference between the two configs.

## How the pipeline works (before this change)

```
submit_edges_sims.py
  → sbatch → pairrunner_edges.py  (one job per workdir)
                → simulator_edges.runconnectedpair_prefire_from_edges()
                    reads:  <workdir>/prefire_simulation_config.json
                              "network" → dhuruva_circuit_config.json   ← circuit config
                    patches (in-memory temp copy): node_sets_file, output_dir
                    writes: <workdir>/bluecellulab_results/rho.h5       ← hardcoded subdir

plot_compare_ndamus_bluecellulab.py
  → process_results(..., rho_subdir="bluecellulab_results")             ← hardcoded in main()
      reads: <workdir>/bluecellulab_results/rho.h5
```

Two blockers:

1. `simulator_edges.py` never patches `"network"` in the temp config → always uses the old circuit.
2. `simulator_edges.py` hardcodes `"bluecellulab_results"` as the output subdir → ion-channels
   run would silently overwrite existing chindemi results.

## Code changes (4 files, all additive)

### 1. `plastyfire/simulator_edges.py`

**`_run_prefire_edges_process()`** — add two new keyword params at the end of the signature:
```python
circuit_config=None,
bcl_subdir="bluecellulab_results",
```

In the temp-config patching block (right after `cfg["node_sets_file"] = node_sets_local`),
add one line:
```python
if circuit_config:
    cfg["network"] = circuit_config
```

Where `bcl_out` is set, replace the hardcoded string `"bluecellulab_results"` with `bcl_subdir`:
```python
# was:
bcl_out = os.path.join(workdir, "bluecellulab_results")
# becomes:
bcl_out = os.path.join(workdir, bcl_subdir)
```
(The `bluecellulab_output_dir` branch already constructs its own path; only the else-branch
uses the hardcoded name.)

**`runconnectedpair_prefire_from_edges()`** — same two new keyword params added to signature.

Apply the same `circuit_config` patch to the snap-read temp config (used only to extract
`pre_gid`/`post_gid`/`t_end` via bluepysnap before the subprocess starts).

Pass both new params to the subprocess via `kwargs=`:
```python
proc = multiprocessing.Process(
    target=_run_prefire_edges_process,
    args=(...existing positional args...),
    kwargs={"circuit_config": circuit_config, "bcl_subdir": bcl_subdir},
)
```

### 2. `plastyfire/pairrunner_edges.py`

Add two optional CLI args:
```
--circuit-config   Path to circuit config JSON (default: use whatever is in prefire_simulation_config.json)
--bcl-subdir       Output subdir name inside workdir for rho.h5 (default: bluecellulab_results)
```

Pass both through to `sim_mod.runconnectedpair_prefire_from_edges()`.

### 3. `submit_edges_sims.py`

Add two optional CLI args:
```
--circuit-config   Forwarded verbatim to pairrunner_edges.py
--bcl-subdir       Forwarded verbatim to pairrunner_edges.py
```

Update `build_runner_args()` to append `--circuit-config` and `--bcl-subdir` when provided.

Auto-default `--bcl-subdir` to `bluecellulab_results_ion_channels` when `--circuit-config`
is provided and `--bcl-subdir` was not explicitly set (same pattern as the `fitted`/`fitting2`
auto-routing).

### 4. `plot_compare_ndamus_bluecellulab.py`

Add one optional CLI arg:
```
--bcl-subdir   Subdir name to look for rho.h5 (default: bluecellulab_results)
```

Pass it to the `process_results()` call in `main()` — that function already accepts
`rho_subdir` as a parameter; only the `main()` call hardcodes it.

## Running the ion-channels variant

```bash
# Step 1: submit simulations (same workdirs, new output subdir)
python submit_edges_sims.py \
    --params chindemi \
    --freq 10Hz \
    --fastforward 280000 \
    --circuit-config /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_ion_channels_circuit_config.json \
    --bcl-subdir bluecellulab_results_ion_channels

# Step 2: generate STDP curves from the new results
python plot_compare_ndamus_bluecellulab.py \
    --bcl-subdir bluecellulab_results_ion_channels \
    --output stdp_compare_10Hz_ion_channels.png
```

## What is NOT touched

| Item | Status |
|---|---|
| All `prefire_simulation_config.json` files on disk | Unchanged — `network` is patched only in the in-memory temp copy |
| Existing `bluecellulab_results/rho.h5` outputs | Untouched — ion-channels run writes to `bluecellulab_results_ion_channels/` |
| Chindemi params / PARAM_PRESETS | Identical |
| Workdir structure / `find_workdirs()` | Same tree, same workdirs reused |
| Default behaviour (no `--circuit-config` flag) | Identical to before this change |
