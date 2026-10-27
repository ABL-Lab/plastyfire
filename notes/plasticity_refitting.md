stdp plot:

```livescript
python plot_stdp_ic_best.py --param-hash 283d79c04294 \
    --label "c_pre-only stride-1 (gamma_p railed)" --out stdp_cpreonly_s1_bcl.png
```

```livescript
cd /lustre06/project/6077694/dhuruva/plastyfire/analytical_method

python fit.py --cpre-only --fit-gamma --de \
              --workers 60 --maxiter 600 --stride 1 --fit-dt=-10,5,10 --anchor-tails

```

# Methods to rerun BCL simulation

## Step 1: Delete old BCL simulation results

These are the `bluecellulab_results_ion_channels/` folders inside every workdir:

```bash
find /project/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations \
    -type d -name "bluecellulab_results_ion_channels" \
    -exec rm -rf {} + 2>/dev/null || true
```

Also delete all the hash-based BCL pkl files:

```bash
find /project/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations \
    -name "simulation_edges_*.pkl" -o -name "simulation_edges.pkl" \
    | xargs rm -f
```

## Step 2: Delete old Neurodamus results (`out_ic/`)

```bash
find /project/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations \
    -type d -name "out_ic" \
    -exec rm -rf {} + 2>/dev/null || true
```

## Step 3: Delete old basis results and recompute

```bash
rm -rf /project/ctb-emuller/dhuruva/plastyfire/basis_results_edges_ion_channels/
```

Then rerun the basis (uses `--circuit-config` pointing to the ion channels config, which picks up the updated hoc):

```bash
cd /project/ctb-emuller/dhuruva/plastyfire

python submit_basis_edges.py \
    --output-dir basis_results_edges_ion_channels \
    --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
    --execution-mode slurm
```

## Step 4: Recompute thresholds from cache (no new BCL needed)

```livescript
cd /project/ctb-emuller/dhuruva/plastyfire

# Rebuild the IC cache against the CURRENT hoc. --circuit-config is required:
# without it the mini-sims read whatever network is baked into each workdir config.
rm -f cpre_cpost_cache/ion_channels_tau278.pkl

python precompute_cpre_cpost.py \
    --params chindemi_aparams \
    --results-dir /project/ctb-emuller/dhuruva/plastyfire/refitting_results \
    --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
    --output cpre_cpost_cache/ion_channels_tau278.pkl \
    --workers 60

```

The thresholds use `compute_thresholds_from_cache.py` which reads the existing cpre/cpost cache — the hoc change doesn't affect c_pre/c_post since those are driven by the GluSynapse dynamics (not the cell ion channels). So you can re-inject with the same a-params:

```bash
cd /project/ctb-emuller/dhuruva/plastyfire

python compute_thresholds_from_cache.py \
    --cache cpre_cpost_cache/ion_channels_tau278.pkl \
    --output-dir threshold_results_ic_gen4/ \
    --a00 1.002     --a01 2.254638 \
    --a10 1.209857  --a11 2.396290 \
    --a20 1.127     --a21 2.500181 \
    --a30 4.214614  --a31 2.239758 \
    --edges-h5 data/dhuruva_modified_edges.h5

python inject_thresholds.py \
    --results-dir threshold_results_ic_gen4/ \
    --edges-h5 data/dhuruva_modified_edges.h5
```

## Step 5: Re-run BCL STDP simulations

```bash
cd /project/ctb-emuller/dhuruva/plastyfire

bash submit_ic_best_700.sh

python run_ic_best_pool.py --workers 60
```

This reads the updated hoc from `modified_emodels_hoc/` (via `dhuruva_modified_ion_channels_circuit_config.json`) and writes to `simulation_edges_0ce64fa83b85.pkl` in each workdir.

```livescript

python plot_stdp_ic_best.py \
    --param-hash cdf3a1e1db98 \
    --basis-dir basis_results_edges_ion_channels \
    --label "IC gen-4, tied apical=basal (new ion channels)" \
    --out stdp_ic_best_tied.png
```

---

**Order matters**: run steps 1–3 first (cleanup + basis), then 4 (thresholds), then 5 (BCL STDP). Steps 3 and 5 are Slurm jobs so you can submit them and wait for completion before the next step. Step 4 runs locally and is fast.

b8c7ff3ecf0a

# Method to create refitting_results folder

with the current set of understanding w.r.t to tau_eff, can you fit STDP to markram data with Natali’s modified params?

answer is no. 

!image.png

even if I leave the tau_eff as dynamic params I cant fit the stdp curves. we need cicr framework to fit them.

the cicr prams arent matching, the mod file changed the cai_CR behavoiur need to fix the mod file first

```python
 python plastyfire/modelfitter_edges.py \
     --basis-dir basis_results_edges_mini/ \
     --cpre-cpost-cache cpre_cpost_cache.pkl \
     --max-workers 62 \
     --pop_size 16 \
     --gen 100 \
     --sample_size 100 \
     -v
```

```livescript
python plastyfire/modelfitter_edges.py \
    --basis-dir basis_results_edges_mini/ \
    --cpre-cpost-cache /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/chindemi_tau278.pkl \
    --max-workers 62 \
    --pop_size 16 \
    --gen 100 \
    --sample_size 100 \
    --seed-individual "101.5,216.2,1.002,1.954,1.159,2.483,1.127,2.456,5.236,1.782" \
    -v

```

For new ion channel 

```python
	python plastyfire/modelfitter_edges.py \
	    --basis-dir basis_results_edges_ion_channels_cached/ \
	    --cpre-cpost-cache /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/ion_channels_tau278.pkl \
	    --circuit-config /project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_ion_channels_circuit_config.json \
	    --max-workers 60 \
	    --pop_size 16 \
	    --gen 100 \
	    --sample_size 30 \
	    --seed-individual "101.5,216.2,1.002,1.954,1.159,2.483,1.127,2.456,5.236,1.782" \
	    -v
	
```

with same params

```livescript
python plastyfire/modelfitter_edges.py \
    --basis-dir basis_results_edges_ion_channels/ \
    --cpre-cpost-cache /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/ion_channels_tau278.pkl \
    --circuit-config /project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_ion_channels_circuit_config.json \
    --max-workers 62 \
    --pop_size 16 \
    --gen 100 \
    --sample_size 30 \
    --seed-individual "101.5,216.2,1.002,1.954,1.159,2.483,1.127,2.456,5.236,1.782" \
    --same-apical-basal \
    -v
```

```livescript
cd /project/ctb-emuller/dhuruva/plastyfire

python plot_edges_fitness.py \
    --checkpoint checkpoint_edges.pkl \
    --out fitness_ic_samebasal.png

```

```livescript
python plastyfire/modelfitter_edges.py \
    --basis-dir basis_results_edges_ion_channels_cached/ \
    --cpre-cpost-cache /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/ion_channels_tau278.pkl \
    --circuit-config /project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_ion_channels_circuit_config.json \
    --max-workers 20 \
    --pop_size 16 \
    --gen 100 \
    --sample_size 30 \
    --seed-individual "101.5,216.2,1.002,1.954,1.159,2.483,1.127,2.456,5.236,1.782" \
    --same-apical-basal \
    -v

```

Cmd to delete and start afresh

```livescript
rm -f /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/checkpoint_edges.pkl* \
      /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache.pkl* \
      /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/logs/tagged_*.csv

find /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting \
    -name "simulation_edges_*.pkl" -delete

rm -rf /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/.cache

```

this is precomputing cpre_cpost cache

```livescript
rm /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/chindemi_tau278.pkl

python precompute_cpre_cpost.py \
    --params chindemi_aparams \
    --results-dir /project/ctb-emuller/dhuruva/plastyfire/refitting_results \
    --output /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/chindemi_tau278.pkl \
    --workers 60
```

thiscmd is for validating cpore cpost:

**What `validate_cpre_cpost.py` now checks:**

- Runs the fixed `_find_cpre_cpost` (IClamp dur=3ms) on 3 test pairs
- For each synapse, computes `theta_d_new = a*c_pre + b*c_post` and `theta_p_new` (using the correct basal/apical a-coefficients)
- Compares both against the `theta_d`/`theta_p` stored in `edges.h5` (which `compute_thresholds.py` built with IClamp)
- Prints per-synapse `td_ratio` and `tp_ratio` — expects 0.8–1.2 for PASS
- Also prints `c_post` values directly, so you can see if apical synapses now have non-zero c_post

**To run after precomputing:**

`sbatch submit_validate_cpre_cpost.sh`

This code is for running BCL prefire simulation inorder to validate with neurodamus

```python
cd /project/ctb-emuller/dhuruva/plastyfire
BCL_SUBDIR_OVERRIDE=bluecellulab_results_edges_chindemi_params_prefire bash submit_bcl_prefire_chindemi.sh

```

For new ion channles:

```livescript
cd /project/ctb-emuller/dhuruva/plastyfire

python submit_basis_edges.py \
    --output-dir basis_results_edges_ion_channels \
    --circuit-config /project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_ion_channels_circuit_config.json

```

```livescript
cd /project/ctb-emuller/dhuruva/plastyfire/optimizer_ion_channels

python ../plot_edges_fitness.py \
    --checkpoint checkpoint_edges.pkl \
    --out fitness_ion_channels.png

```

```livescript
python plot_stdp_ic_chindemi.py
```

```livescript
python plot_stdp_ic_best.py

```

## Commands to rewrite `dees_edges.h5` for cell 186168

Run these from `/home/dhuruva/projects/ctb-emuller/dhuruva/DEES_cell_packages/cpre_cpost` on a compute node (per your usual practice, not the login node):

```bash
cd /home/dhuruva/projects/ctb-emuller/dhuruva/DEES_cell_packages/cpre_cpost

# 1. Clear the stale CPre cache — run_cpre_parallel.py skips synapses already in this CSV
rm -f dees_cpre_values_csv_DEES_gmax_rho_50_50/dees_cpre_values_186168.csv

# 2. Recompute CPre for all synapses (writes to dees_cpre_values_csv_DEES_gmax_rho_50_50/)
python run_cpre_parallel_specific_synapses.py     --post_cell_id 186168     --synapse_ids_csv scratch_selections/selected_synapses_186168_loc_35_43_44.csv     --num_workers 60 --batch_size 5

# 3. Recompute CPost (overwrites unconditionally, no cache guard)
python compute_dees_cpost.py --post_cell_id 186168

# 3. Backup before in-place edit (recommended)
cp ../DEES_gmax_rho_50_50/186168/dees_edges.h5 \
   ../DEES_gmax_rho_50_50/186168/dees_edges.h5.bak_pre_ionchannel_update

# 4. Write theta_d/theta_p into dees_edges.h5 — updates every synapse present in
#    BOTH the cpre and cpost CSVs, i.e. only your 83 targeted synapses this time,
#    since dees_cpre_values_186168.csv will only contain those 83
python change_thresholds.py
```

Notes:

- Step 2's `-num_workers`/`-batch_size` are tunable; 18 workers assumes you have that many cores on the node — adjust to what you request from SLURM.
- `change_thresholds.py` has `post_cell_id = 186168` hardcoded and no CLI args — it's already scoped correctly, just run it as-is.
- It edits `DEES_gmax_rho_50_50/186168/dees_edges.h5` **in place** (`r+` mode) — make a backup first if you want to be able to diff/revert:
    
    ```bash
    cp /home/dhuruva/projects/ctb-emuller/dhuruva/DEES_cell_packages/DEES_gmax_rho_50_50/186168/dees_edges.h5 \
       /home/dhuruva/projects/ctb-emuller/dhuruva/DEES_cell_packages/DEES_gmax_rho_50_50/186168/dees_edges.h5.bak_pre_ionchannel_update
    ```
    
- After this, `submit_job_loc.py` will pick up the new thresholds automatically since it reads `dees_edges.h5` live via `SSim`/`instantiate_gids`.

# Running plasticity sims Neurodamus:

Here's what was set up and how to run the full validation:

**What was changed:**

- `submit_edges_sims.py` — added `ic_gen4` preset (gamma_d=101.5, gamma_p=199.773931, tau=278.32)
- `submit_neurodamus_ic.py` — new script, 700 workdirs already have `prefire_simulation_config_ic.json` + `neurodamus_sbatch_ic.sh` generated
- `plot_compare_ndamus_bluecellulab.py` — added `-nd-subdir` arg, fixed `sort_elements` to handle `out_ic`

**To run the validation:**

```bash
# 1. Neurodamus side — 700 jobs writing out_ic/rho.h5
python submit_neurodamus_ic.py

# 2. Bluecellulab side — 700 jobs writing bluecellulab_results_ion_channels/rho.h5
python submit_edges_sims.py \
    --params ic_gen4 \
    --freq 10Hz \
    --fastforward 280000 \
    --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
    --cpre-cpost-cache cpre_cpost_cache/ion_channels_tau278.pkl \
    --skip-existing

# 3. Compare once jobs finish
python plot_compare_ndamus_bluecellulab.py \
    --nd-subdir out_ic \
    --bcl-subdir bluecellulab_results_ion_channels \
    --basis-dir basis_results_edges_ion_channels \
    --output stdp_compare_ic_gen4.png
```

All 700 jobs are submitted — canary 64842543 running, the rest queued (Priority). Each has a 2-hour limit. Results will land in `out_ic/rho.h5` per workdir.

Once the canary finishes and confirms the sim completes, you'll have all 700 running in parallel. After they're done, the next step is the bluecellulab side:

```bash
python submit_edges_sims.py --params ic_gen4 --freq 10Hz --fastforward 280000 \
  --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
  --cpre-cpost-cache cpre_cpost_cache/ion_channels_tau278.pkl --skip-existing
```

Then once both are done:

```bash
python plot_compare_ndamus_bluecellulab.py --nd-subdir out_ic \
  --bcl-subdir bluecellulab_results_ion_channels \
  --basis-dir basis_results_edges_ion_channels \
  --output stdp_compare_ic_gen4.png
```