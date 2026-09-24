# Reproducing the STDP fit (Aug 2026)

Fit the Graupner-Brunel a-params to Markram 1997 in-vitro STDP for L5TTPC→L5TTPC
pairs, using an offline analytical model instead of a simulation-in-the-loop GA.

Repo root `$PF = /project/ctb-emuller/dhuruva/plastyfire`. **Run everything on a
compute node**, never the login node.

**Result (DE fit #2 — current best):**
`a00=1.003498 a01=2.902478 a10=1.644558 a11=2.764812` (apical tied to basal),
`gamma_d=77.7558 gamma_p=299.9121`, `tau_effca=278.3177658387` → hash
`b8c7ff3ecf0a`. Both gammas are inside the physiological bounds
`gamma_d ∈ [50, 200]`, `gamma_p ∈ [150, 300]`.

| Δt | in vitro | bluecellulab (658 pairs) |
|---|---|---|
| −10 | 0.7922 | **0.8065** |
| +5 | 1.2038 | **1.2376** |
| +10 | 1.2013 | **1.1549** |

Full BCL sweep: −50 0.8435, −30 0.8423, −10 0.8065, +5 1.2376, +10 1.1549,
+30 0.9570, +50 0.8403.

The offline model's own prediction for this parameter set was not written to a
log, so only the measured column is recorded.

**Superseded — DE fit #1:** `a00=1.100479 a01=2.891376 a10=1.543869 a11=3.095199`,
`gamma_d=139.6809 gamma_p=483.4465` → hash `9e43ab966326`
(BCL: −10 0.8287, +5 1.1703, +10 1.0954; offline: 0.8112 / 1.1953 / 1.1791).
Discarded because `gamma_p=483` is outside the physiological range. Its
neurodamus jobs were cancelled mid-flight; 379 of 700 `out_defit/rho.h5` from
that run remain on disk and are **not** comparable to `out_defit2/`.

---

## 0. Prerequisites — code state

These are required; without them the numbers below do not reproduce.

1. **`configs/L5TTPC_L5TTPC_STDP.yaml:6`** → `data/dhuruva_modified_ion_channels_circuit_config.json`
   (was `dhuruva_circuit_config.json`, i.e. the unmodified emodels).
2. **`plastyfire/simwriter.py::check_electrical_constraint`** — searches a
   single-AP threshold in addition to the 5-spike induction one, stores it under
   key `"single_ap"`, and **rejects the cell** if none is found. Without the
   rejection, cells with no valid c_post enter the pair list and fail hours later.
3. **`plastyfire/simulator_edges.py`** — c_post uses `spike_threshold_finder`
   with widths `[1.5, 3, 5]` ms over `[0.05, 5.0]` nA / 100 levels (matching
   `thresholdfinder.py`); plus an induction guardrail that raises if the post
   cell does not fire exactly `nspikes × nreps = 50` spikes.

**The mechanism actually loaded is NOT `plastyfire/GluSynapse.mod`.**
`MECHANISMS_PATH` points at `DEES_cell_packages/`, whose compiled
`x86_64/GluSynapse.cpp` came from `other_mods/modified_mechanisms/GluSynapse.mod`
and contains

```
- dep_GB*(1 - pot_GB)*gamma_d_GB*rho_GB
```

The `(1 - pot_GB)` factor makes depression mutually exclusive with potentiation
(a real depression window). The repo's `.mod` lacks it. Any offline model must
include it.

---

## 1. Simulation workdirs

```bash
cd $PF
mv refitting_results refitting_results.old_$(date +%Y%m%d)      # if one exists
python -c "
from plastyfire.simwriter import OptSimWriter
w = OptSimWriter('configs/L5TTPC_L5TTPC_STDP.yaml')
w.write_sim_files(w.find_pairs())"
```

Must run from `$PF` (`templates/simulation.batch.tmpl` is a bare relative path).
30-way parallel NEURON.

→ 100 pairs × 7 Δt = **700 workdirs**, 107 `single_cells/*.pkl` (105 with a
`single_ap` entry, all found at width 1.5 ms, amp 0.80–2.35 nA).

Verify:
```bash
python -c "
import glob,pickle
fs=glob.glob('refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations/single_cells/*.pkl')
print(sum('single_ap' in pickle.load(open(f,'rb')) for f in fs), '/', len(fs))"
```

## 2. EPSP basis

```bash
python submit_basis_edges.py \
    --output-dir basis_results_edges_ion_channels \
    --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
    --execution-mode slurm --walltime 00:30:00 --mem 10g --skip-existing
```

100 SLURM jobs, 12 cores each. Measured: 8–22 min, MaxRSS 7.8 GB at `--mem 10g`
(tight — 22% headroom). → **100 `basis_*.csv`**.

## 3. c_pre / c_post cache

```bash
rm -f cpre_cpost_cache/ion_channels_tau278.pkl
python precompute_cpre_cpost.py \
    --params chindemi_aparams \
    --results-dir $PF/refitting_results \
    --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
    --output cpre_cpost_cache/ion_channels_tau278.pkl --workers 60
```

`chindemi_aparams` is correct despite the name: only `tau_effca` (278.3177658387)
shapes `effcai`. `gamma_d`/`gamma_p` live in the rho ODE, downstream of calcium,
and never touch the cached values.

→ **98/100 pairs, 714 synapses**, ~25 min. Two pairs fail (`189500-185042`,
`194685-197123`) — no single-AP stimulus, so no valid c_post. 22.1% of synapses
have c_post at the resting floor (<1e-4).

Steps 2 and 3 are independent and can run concurrently.

---

# DE fitting (analytical method)

Bootstraps off one existing bluecellulab run: it needs `simulation_traces.pkl`
(calcium) and `simulation_edges_*.pkl` (rho ground truth). Run the **BCL
simulation** section once with any parameter set first — the tied-ic4 run
(`cdf3a1e1db98`) was used here.

## 4. Extract calcium traces

```bash
cd $PF/analytical_method
python extract.py --workers 60
```

Reads 152 GB of traces, rebuilds `effcai` at native 0.025 ms
(`effcai' = -effcai/tau_effca + (cai_CR - min_ca_CR)`, exact exponential
integrator), decimates to 2 ms.

→ **653 `.npz` in `extracted/`, 381 MB**. Nothing here depends on the a-params —
run once, then every candidate is evaluated from these arrays.

*If the cache is lost later, recover it from these without re-simulating:*
```bash
python rebuild_cache.py          # exact; valid only while tau_effca is unchanged
```

## 5. Validate the model against bluecellulab

```bash
python validate.py
```

Gate before fitting anything. Expect:

```
corr(pred, obs)   1.0000     mean |pred-obs|  0.0011
BINARY agreement  99.9%      curve error      max 0.0104
```

3 of 7 Δt reproduce to 4 decimals. The residual is 6 synapses of 4729 sitting
within 0.0035 of the `rho ≥ 0.5` binarisation cut (`python residual.py` confirms).
Model error is ~7× below the measurement SEM.

## 6. Fit

```bash
python fit.py --time-eval --fit-gamma --workers 60        # ~2.9 s/eval
python fit.py --fit-gamma --de --workers 60 --maxiter 80 \
              --anchor-tails --sign --verbose
```

~3–5 min on 60 cores. Objective = Σ((pred−target)/SEM)² at Δt = −10, +5, +10.

**Why DE/CMA-ES and not gradients or a NN:** the EPSP stage binarises rho at 0.5,
so the objective is piecewise-constant — flat almost everywhere, jumping only
when a synapse crosses the cut (4682 of 4729 synapses are far from it and
contribute no local sensitivity). No gradient exists. A surrogate is pointless
when an evaluation costs 3 s. `--de` is scipy differential evolution; drop it for
CMA-ES (`pip install --user cma`) and treat agreement between the two as evidence.

Four constraints matter — without them the search finds degenerate optima:

| rule | why |
|---|---|
| `θ_p > θ_d` per synapse | coefficient tests (`a10>a00`) are neither necessary nor sufficient; θ depends on each synapse's own `(c_pre, c_post)` |
| ≥30% synapses with `θ_d < peak effcai` | else nothing crosses, rho frozen, every Δt = 1.0 (a huge flat basin) |
| ≥15% with `θ_p < peak effcai` | else `pot` never fires → "depress everything": fits −10 while killing +5/+10 |
| a-params ∈ [1, 5] | unbounded search ran `a01→37`, `a11→47`, putting θ at 2× the calcium |

`--anchor-tails` pulls ±30/±50 toward 1.0 (loose, SEM 0.15 — not measured points).
`--sign` penalises wrong-sign plasticity outright.

→ weighted err **2.594** in 3888 evaluations (from 544.67 at the seed).

Parallelism note: `Batch` is a module-level global built before the pool forks, so
60 workers share ~800 MB copy-on-write. The objective must stay a top-level
function — as a closure it gets pickled and the sharing is lost.

---

# BCL simulation (bluecellulab)

θ is computed **at runtime** from the cache + a-params
(`use_a_params` branch in `simulator_edges`), so `inject_thresholds.py` is **not
required**. `edges.h5` `theta_d`/`theta_p` are ignored on this path — verified:
686 log lines "Loaded c_pre/c_post from cache", and runtime θ matched the injected
values on 714/714 synapses.

```bash
cd $PF
python run_de_fit2_pool.py --workers 60          # DE fit #2,  hash b8c7ff3ecf0a  <- current best
# or: run_de_fit_pool.py                          # DE fit #1,  hash 9e43ab966326  (superseded)
# or: run_ic_best_pool_tied.py                    # tied ic4,   hash cdf3a1e1db98
```

~4 GB per worker (60 workers ≈ 240 GB). 700 tasks, ~20 min.
→ `simulation_edges_<hash>.pkl` per workdir, plus `bluecellulab_results_optimizer/simulation_traces.pkl`
(the calcium traces step 4 consumes).

Plot:
```bash
python plot_stdp_ic_best.py --param-hash b8c7ff3ecf0a \
    --label "DE fit2 (analytical)" --out stdp_de_fit2.png
```

**`--param-hash` is mandatory** — the script defaults to `0ce64fa83b85` and
silently collects zero ratios otherwise.

Expect ~658/700 pkls. The shortfall is 2 uncached pairs × 7 = 14, plus ~33
induction-guardrail failures (51–60 post spikes instead of 50, clustered at short
|Δt| where the EPSP summates with the pulse). Every completed run delivers
**exactly 50/50** post spikes.

### Adding a new parameter set

Copy `run_de_fit_pool.py`, change `FIT_ARGS` **and** `PARAM_HASH`. The hash names
the output file, so changing params without the hash silently mixes two runs.

> **`b8c7ff3ecf0a` is a non-canonical label.** It is the md5 of the discarded
> `gamma_p=738.83` parameter set; `FIT_ARGS` was later refit under the
> physiological bounds without updating `PARAM_HASH`. The canonical hash of the
> params actually in `run_de_fit2_pool.py` is `1f1909bc65dc`. No data was mixed —
> the 738.83 set was never simulated, and all 658 `simulation_edges_b8c7ff3ecf0a.pkl`
> post-date the refit — but the hash no longer identifies its params. Fit #1's
> `9e43ab966326` does verify against the snippet below.

```bash
python -c "
import hashlib
fp={'gamma_d_GB_GluSynapse':77.7558,'gamma_p_GB_GluSynapse':299.9121,
    'a00':1.003498,'a01':2.902478,'a10':1.644558,'a11':2.764812,
    'a20':1.003498,'a21':2.902478,'a30':1.644558,'a31':2.764812,
    'tau_effca_GB_GluSynapse':278.3177658387}
print(hashlib.md5(str(sorted(fp.items())).encode()).hexdigest()[:12])"
```

---

# Neurodamus

Different from BCL in two ways. Both steps are **mandatory**.

Neurodamus has no runtime a-param path — it reads `theta_d`/`theta_p` straight
from `dhuruva_modified_edges.h5`, and takes the gammas from
`conditions.mechanisms.GluSynapse` in the config.

## 1. Inject the thresholds

```bash
cd $PF
python compute_thresholds_from_cache.py \
    --cache cpre_cpost_cache/ion_channels_tau278.pkl \
    --output-dir threshold_results_de_fit2/ \
    --a00 1.003498 --a01 2.902478 --a10 1.644558 --a11 2.764812 \
    --a20 1.003498 --a21 2.902478 --a30 1.644558 --a31 2.764812 \
    --edges-h5 data/dhuruva_modified_edges.h5

python inject_thresholds.py \
    --results-dir threshold_results_de_fit2/ \
    --edges-h5 data/dhuruva_modified_edges.h5 --force
```

`--force` is required here: `edges.h5` still holds the DE fit #1 thresholds, and
`inject_thresholds.py` aborts on a mismatch unless told to overwrite.

Rewrites the shared 4.2 GB `edges.h5` **in place** and is not versioned. Verify
714/714 synapses carry `theta_p > theta_d`.

## 2. Run

```bash
python submit_neurodamus_defit.py --dry-run
python submit_neurodamus_defit.py
```

Sets `gamma_d=77.7558`, `gamma_p=299.9121` in the config; writes
`prefire_simulation_config_defit2.json`, `out_defit2/`, jobs `nddefit2_*` — no
collision with `out_ic/` or with DE fit #1's `out_defit/`. All four names derive
from the single `TAG` constant at the top of the script. A startup guard recomputes θ from the cache and aborts
if `edges.h5` disagrees on >0.1% of synapses (`--skip-threshold-check` overrides).

---

## Known caveats

- **`gamma_p` sits at its ceiling** (299.91 of a 300 bound). The search wanted
  more potentiation than the physiological range allows — first unbounded
  attempts railed at 483 then 739. Read any future demand for higher `gamma_p` as
  the model under-potentiating, not as a better fit.
- **+10 ms is the weakest point**: 1.1549 measured vs 1.2013 in vitro. DE fit #1
  showed the same sign of error (1.0954), and there the offline/BCL gap was 0.084
  against a 0.0104 validation error.
- **Tails stay depressed**: −50 0.8435, −30 0.8423, +30 0.9570, +50 0.8403, where
  physiology expects ≈1.0. Most of the residual fit error is exactly this. Not
  yet explained.
- **Apical/basal is tied** (`a20=a00` etc). An untied 8-param fit needs
  per-synapse section labels (`sec_types[gid] == 3`), which `extract.py` does not
  save.
- **22.1% of synapses have c_post at the resting floor**, so their depression
  window `(a11−a01)·c_post ≈ 0` regardless of parameters. A structural ceiling on
  achievable LTD.
- `plastyfire/GluSynapse.mod` is **not** the running mechanism (see §0).
