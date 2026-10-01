# kimchi for the GluSynapse_v2 GPU fits (plan, 2026-09-30)

Research and plan only. Nothing has been pushed, created on GitHub or submitted through kimchi. Draft repo: `/lustre09/project/6070394/dhuruva/glusynapse-fits_draft/` (copies only).

## 1. What kimchi is

Source: `/lustre09/project/6070394/dhuruva/kimchi` (README.md, `kimchi/*.py`, `kimchi/job.sh`). Operational notes: `kettle/CLUSTER_RUNBOOK.md`.

- **Where it runs.** kimchi runs on the **nibi login node** only. Binary: `~/.local/bin/kimchi` → `/project/6055588/dhuruva/kimchi/.venv/bin/kimchi`. Ledger: `/project/6055588/dhuruva/kimchi_ledger`, holding `kimchi.toml`, `events.jsonl` and `results/`.
  - It is **not installed on rorqual**: `which kimchi` finds nothing, and there is no `~/.config/kimchi` and no `~/.ssh/config`. The copy here is the source checkout plus rorqual's spool (`/lustre09/project/6070394/dhuruva/kimchi_spool`).
  - I could not read the real `kimchi.toml` from rorqual. The accounts and flags below come from the docs, marked as such.
- **How it works.** kimchi races one run on every eligible cluster. The first copy to *start* claims the run, and kimchi cancels the other copies:
  - pending copies are always cancelled;
  - running copies are cancelled only while younger than `grace_seconds`, default 600 s.
- **Remote transport.** Remote clusters are reached over ssh ControlMaster sessions that `kimchi login` opens: one Duo push per cluster, then `BatchMode`.
  - Nothing is installed remotely. `agent.py` (stdlib only) is piped over ssh.
  - A cluster whose session is dead is marked `UNKNOWN`. kimchi never submits to it or cancels anything on it.
- **Code transport is git.** `kimchi submit` does `git push` from nibi's checkout, then `git fetch` + `--ff-only` in each cluster's checkout (`[projects.X] dirs`).
  - It races only on clusters that end up on the same commit with a clean tree.
  - **Data is not moved by kimchi.** It must already be on each cluster, which is why we use Globus.
- **Sizing.** Per cluster, `kimchi.toml` maps `gpus = { <type> = "<sbatch flags for ONE gpu>" }`. `{bucket}` is filled with the b1..b5 walltime bucket.
  - A sweep chooses `gpu = [...]` (a preference list), `time`, and `mem`. `mem` is appended as `--mem=` and overrides the cluster flag.
  - **CPUs cannot be set per sweep.** So our project needs its own GPU key, with `--cpus-per-task=1` and a MIG slice (see 3).
  - Whole-node `pack` clusters (tamia) are opt-in, and they pack N runs per node.
- **Results.** The job wrapper (`job.sh`) writes started/finished records to `<spool>/inbox`. On `sync`, kimchi reads the plugin's `results_path` JSON over ssh and stores it in `<ledger>/results/<sweep>/<run>.json`. That file carries provenance: cluster, job, node, commit, and a frozen copy of the sweep.
  - Only the JSON (the `KEEP` keys) comes back. Other outputs (csv, ckpt, logs) stay on the cluster.
  - Logs are in `<spool>/logs/` on the winning cluster.
- **Waiting and status.** Nothing notifies you: the race advances only while something polls.
  - `kimchi submit <sweep> --watch` (or `kimchi watch`, in tmux on nibi) polls every 60 s while copies are queued, then every `poll_seconds` (default 300).
  - `kimchi status <sweep> --runs` does one refresh and returns.
  - `kimchi status --live` only re-reads the ledger (no ssh).
  - Agents should run a single `kimchi status` and must not loop on it.

### Commands

| command | effect |
|---|---|
| `kimchi login [clusters…]` | opens the ssh sessions, one Duo prompt per cluster; `--keepalive` reconnects |
| `kimchi show <sweep>` | shows the resolved sweep and its run count |
| `kimchi submit <sweep> [--dry-run] [--watch] [--project glusyn]` | pushes code, fast-forwards the clusters, sbatches the missing runs. **The only command that starts work.** Always name the sweep: the ledger is shared with MMCR_CL/kettle |
| `kimchi status [sweep] [--runs] [--offline] [--live] [--json]` | refreshes, then shows pending/running/done per cluster |
| `kimchi sync` | one refresh: collects results and cancels losers; starts nothing |
| `kimchi results <sweep> [--json]` | metric table (our metric is chi2) |
| `kimchi stop <sweep> [--undo]` / `kimchi retry <sweep> <run>` | cancels pending copies / re-races one run |
| `kimchi adopt <sweep> jobs.tsv` | tracks jobs that were submitted by hand |

## 2. Clusters and GPUs

Sources: `kettle/CLUSTER_RUNBOOK.md` and `MMCR_CL/mmcr_env.py`. The accounts are MMCR's defaults; the real `kimchi.toml` on nibi is authoritative.

| cluster | GPU kimchi uses | MIG? | account (docs) | kettle/MMCR checkout root |
|---|---|---|---|---|
| nibi (kimchi host) | h100 | verify (`sinfo -o %G`) | def-emuller_gpu | /project/6055588/dhuruva |
| fir | h100 | verify | def-emuller_gpu | /scratch/dhuruva |
| rorqual (here) | h100 | **yes**: 1g.10gb / 2g.20gb / 3g.40gb (our pilots ran on them) | def-emuller (→ _gpu) | /project/rrg-emuller/dhuruva |
| vulcan | l40s 48 GB | no (L40S has no MIG) | aip-emuller | /project/aip-emuller/dhuruva |
| killarney | h100 or l40s | verify | aip-emuller | /project/6102315/dhuruva |
| tamia | 4× h100 whole node (`pack`, opt-in, min_time 1 h) | no | aip-emuller | /project/6102315/dhuruva |
| trillium / trillium-gpu | login only, `account="TBD"` | — | — | do not use |

## 3. Our fit through kimchi

Measured on rorqual. Pilot 22098228 ran the 5 -vca dirs, paired_l5+sjostrom07 and CHUNK 4 on a 1g.10gb slice:

- MaxRSS 35.9 GB, CPU efficiency 72% on 1 CPU;
- objective call 13–15 s, which extrapolates to about 135 min for MAXITER 150.

That gives the request **1 CPU, `--mem=45G`, `--time=03:30:00`, 1g.10gb MIG**. The four real fits, 22098444–47, are running now. Run `seff` on them when they finish and update `sweeps/_base.toml`.

**kimchi.toml additions** (on nibi; user to confirm the accounts and paths):
```toml
[clusters.rorqual]   # add a key; keep the existing h100 entry for MMCR/kettle
gpus = { h100_1g = "--gpus-per-node=nvidia_h100_80gb_hbm3_1g.10gb:1 --cpus-per-task=1 --mem=45G" }
# same key on nibi/fir only after `sinfo -o %G` shows the MIG gres there.
# vulcan/killarney l40s: a whole 48 GB GPU per fit. Allow only if the user accepts that (see Q3).

[projects.glusyn]
plugin = "kimchi.plugins.glusyn"        # ln -s <checkout>/kimchi_plugin.py <kimchi>/kimchi/plugins/glusyn.py
dirs = { rorqual = "/project/rrg-emuller/dhuruva/glusynapse-fits", nibi = "/project/6055588/dhuruva/glusynapse-fits" }
sweeps = "sweeps"
```

**Worked example.** This is `run_fits_td4_pl5sj07.sh` as the sweep `sweeps/td4_pl5sj07.toml`: td4 seeds 1 and 2 seeded from td3_s5, td4 seed 3 unseeded (the basin check), and the td3 s1 reference. That is 4 runs.
```sh
# on the nibi login node, inside tmux
kimchi login rorqual fir                      # Duo pushes
kimchi show td4_pl5sj07                       # expect 4 runs
kimchi submit smoke --dry-run && kimchi submit smoke    # 10 min / 2G jax.devices() check (sweeps/smoke.toml)
kimchi status smoke --runs                    # one check; the job log is <spool>/logs on that cluster
kimchi submit td4_pl5sj07 --dry-run           # prints the sbatch argv per cluster
kimchi submit td4_pl5sj07 --watch             # races until settled; Ctrl-C is safe
kimchi results td4_pl5sj07                    # chi2 per run
```
- **Each run.** It executes `scripts/run_fit.sh TD=4 SEED=1 SEEDFITS=... ...` and writes `results/<run>/fit.{json,csv}` plus `fit_ckpt.npz`.
- **Resume.** A retry on the same cluster resumes automatically from the ckpt.
- **Login-node safety.** The plugin's `check` is stdlib only and does no subprocess python, so nothing heavy runs on the login node.

## 4. Slim repo proposal: `glusynapse-fits`

The layout **mirrors plastyfire's relative paths**. That way `fit_v2.py`, `batch_v2.py` and `targets.py` run with **no code edits** (they resolve `../analytical_method`, `../ebner` and `../nevian`), and `diff -r` against plastyfire shows any drift.

```
glusynapse-fits/
  glusynapse_v2/  fit_v2.py jax_v2.py batch_v2.py model_v2.py targets.py   (import closure, unchanged copies)
                  edge_params.npz (47 KB, edge cache)  subset24_pairs.txt
  analytical_method/constants.py      # the ONLY file needed from analytical_method (31 GB dir), no imports
  ebner/ebner_targets.csv             # targets.py: paired_l5 silently drops the Sjostrom csv targets without it
  nevian/nevian_protocols.py          # only for the 'nevian' group (stdlib only)
  fit_results/delta-cooker.json       # --x0
  seeds/reduced_gpu_subset_pl5r50_td3_s5.json   # --seed-fits
  scripts/run_fit.sh                  # kimchi job: env-driven port of run_reduced_gpu.sh (MODE=subset)
  scripts/sync_from_plastyfire.sh     # plastyfire -> repo (cmp + cp, never writes plastyfire)
  scripts/results_to_plastyfire.sh    # results/<run>/fit.* -> plastyfire/glusynapse_v2/results/reduced_gpu_subset_<run>.*
  sweeps/_base.toml  td4_pl5sj07.toml  smoke.toml
  kimchi_plugin.py   env/setup_cluster.sh   requirements.txt   .gitignore (results/, data/, *_ckpt.npz)
```

- **Import closure** (found by grep):
  - `fit_v2` → `batch_v2`, `targets`, `model_v2`, `jax_v2`
  - `batch_v2` → `analytical_method/constants`, `model_v2`, numba, pandas
  - `targets` → `nevian_protocols` (lazy), `ebner_targets.csv`
  - `model_v2` → scipy.signal; h5py only on an `edge_params` cache miss
- **Pins** (from `plastyfire/.venv`, Python 3.11.4): numpy==1.26.4, scipy==1.17.1, pandas==2.3.3, numba==0.67.0, llvmlite==0.49.0, jax[cuda12]==0.7.1 (jaxlib/cuda12 plugin 0.7.1), h5py==3.16.0.
  - `env/setup_cluster.sh` builds `~/gsf_venv` on each cluster's login node (install only).
- **Data location.** Data dirs are set by env vars in a per-cluster `~/.gsf_env`, which is not in the repo:
  - `GSF_DATA` is the Globus destination;
  - `ANALYTICAL_BASIS_DIR` defaults to `$GSF_DATA/basis_results_edges_sabrina_n120_delta`;
  - `PLASTYFIRE_EDGES=/nonexistent`, so a cache miss fails loudly;
  - `GSF_VENV`.
  - Sweeps name dirs relative to `$GSF_DATA/extracted`. There are no `/lustre` paths.
- **Two silent-failure traps** to check in the smoke/first run:
  - `batch_v2` **skips pairs with no basis csv without an error**. Check the "skipped" count in the log.
  - `edge_params.npz` must cover every synapse in the npz files, because the h5 is not on the remote clusters. The rorqual pilot ran from this cache, so it should.

## 5. Globus transfer list

Source: rorqual `/project/rrg-emuller/dhuruva/plastyfire/`. Destination: `$GSF_DATA/` on each cluster.

| path (relative to plastyfire) | size | files |
|---|---|---|
| glusynapse_v2/extracted/ebner_delta-prefire-vca | 1.6 GB | 240 npz |
| glusynapse_v2/extracted/markram_delta-prefire-vca | 1.2 GB | 168 |
| glusynapse_v2/extracted/sj03_delta-prefire-vca | 2.9 GB | 220 |
| glusynapse_v2/extracted/sj03r50_delta-prefire-vca | 838 MB | 48 |
| glusynapse_v2/extracted/sj07_delta-prefire-vca | 926 MB | 94 |
| basis_results_edges_sabrina_n120_delta → `$GSF_DATA/basis_results_edges_sabrina_n120_delta` | 7.5 MB | 122 csv |
| **total** | **≈ 7.5 GB** | |

These go through git, not Globus: `fit_results/delta-cooker.json` (13 KB), the seed fit json, `edge_params.npz` (47 KB), `ebner_targets.csv` (13 KB), `subset24_pairs.txt` (13 KB). All of `glusynapse_v2/results/` is 13 MB, but only the named seed fits are needed.

Not needed remotely: the circuit, the simulation outputs, `dhuruva_modified_edges.h5`, `analytical_method/` beyond `constants.py`.

## 6. Results back and staying in sync

- **Small outputs.** The fit json (a, pre, chi2, args, minutes) comes back automatically into the nibi ledger `results/glusyn...`, with provenance.
- **Large outputs.** Use a Globus pull of `<checkout>/results/<run>/` (fit.csv, fit_ckpt.npz) plus the `<spool>/logs` of the winning job. Then, on rorqual, `scripts/results_to_plastyfire.sh <dir>` copies them in under the existing `reduced_gpu_subset_<run>.*` names and never overwrites.
- **Code direction.** plastyfire stays the source of truth for model code. `scripts/sync_from_plastyfire.sh` updates the repo copies, then commit. New fit ideas are developed in the slim repo and ported back to plastyfire by hand once they are kept.
- **Records.** Log each kimchi sweep, with the ledger path and its `seff`, in `glusynapse_v2/DECISIONS.md`.

## 7. Open questions for the user

1. **Repo.** Name and owner: `dhuruvapriyan/glusynapse-fits` (like kettle) or `ABL-Lab/...`? I need your confirmation before creating or pushing anything.
2. **Ledger.** Share the nibi ledger with MMCR_CL/kettle, or `kimchi init` a separate one? A bare `kimchi watch` on the shared ledger also drives their sweeps.
3. **L40S.** On vulcan/killarney, allow a whole L40S per fit (48 GB, no MIG) or restrict to H100 MIG lanes? A whole GPU for a job that fits in a 1g slice (1/7 of an H100) is billed at several times the fairshare.
4. **MIG elsewhere.** Do nibi/fir/killarney offer the `nvidia_h100_80gb_hbm3_1g.10gb` gres? Check with `sinfo -o %G` on each before adding the `h100_1g` key there.
5. **tamia.** Use it for packs of 4 seeds per node? One node = 4× H100, `--mem` = 4 × 45G. Or leave it out.
6. **Accounts.** Confirm the GPU accounts per cluster in the real `kimchi.toml`, and that aip-emuller use is OK for this project.
7. **Globus paths.** Destination endpoints and paths per cluster; these set `GSF_DATA` in `~/.gsf_env`.
8. **Smoke test.** Not run. kimchi is only usable from nibi with live Duo sessions, which this agent does not have. Run `kimchi submit smoke` once the repo is cloned on rorqual.
