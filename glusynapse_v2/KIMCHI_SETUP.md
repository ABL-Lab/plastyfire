# kimchi setup for the stage-4 GPU screens (plan, 2026-10-02)

This is a research and planning document. Nothing has been changed or submitted. It builds on `KIMCHI.md` (2026-09-30), checks that file against the kimchi source and corrects it. The target is about 64 rule variants × 4 starts, so about 256 fit+val runs, raced across clusters.

Sources read: kimchi at `/project/rrg-emuller/dhuruva/kimchi` (README, `kimchi/*.py`, `job.sh`, `plugins/`), `kettle/CLUSTER_RUNBOOK.md`, the glusyn ledger `kimchi_ledger_glusyn/`, both spools, the repo `/lustre09/project/6070394/dhuruva/glusynapse-fits`, and the stage-4 entry points with their import closure.

## Key findings

1. **The stage-4 fit is a different program from the one KIMCHI.md and the repo were built for.** Stage 4 runs `rho_redesign/fit_launch.py → gpu_v10_rho.py / gpu_v11_rho.py → fit_v7`, which uses numba-cuda kernels and never imports jax. It needs `.venv_v3` (numba 0.65.1 + numba-cuda 0.30.2 + pip CUDA 12.9), a 3g.40gb slice, 29–92 GB of host RAM and the split2 data in `/scratch/dhuruva/s2g0321` (about 22 GB). The repo holds the older fit_v2/jax td4 fit: a 1g.10gb slice, 45G and 7.5 GB of `-prefire-vca` data. Of that, only the env scripts, `requirements_v3.txt` and the plugin skeleton are still useful.
2. **kimchi already runs on rorqual, and through its own ledger.** `~/.config/kimchi/home` points to `/lustre09/project/6070394/dhuruva/kimchi_ledger_glusyn`, and rorqual is the local cluster there. A smoke run (job 22105657, `h100_1g`, 2026-09-30) finished and returned its JSON. This ledger is not shared with MMCR_CL or kettle; their ledger is on nibi. So a bare `kimchi watch` on rorqual only drives glusyn sweeps.
3. **The GitHub repo exists, but local `main` is one commit ahead.** The repo `dhuruvapriyan/glusynapse-fits` exists (remote HEAD `8555fbb`). Local `main` at `b6426e1` (bootstrap script, requirements_v3, GLOBUS.md) was never pushed. rorqual's key can reach GitHub without a prompt. `glusynapse-fits_draft/` does not exist: the draft became `glusynapse-fits`.
4. **ssh sessions are live for three clusters only.** ControlMaster sessions for fir, killarney and vulcan have been alive on login node **rorqual3** since 2026-09-30 17:25 (`ps`). nibi and tamia have none. kimchi must be driven from rorqual3, or the user must run `kimchi login` again from whichever login node is used.
5. **Only the JSON comes back, so the run must summarise itself.** Every remote output except one JSON stays on the winning cluster. To make the selection rule mechanical, the run script has to:
   - run val in the same allocation;
   - write one summary JSON with the Markram ±10 rows, χ² per pathway, params, de_fun, the val rows, peak RSS and the GPU name;
   - return the whole fit JSON inside it, so that val or LOPO can be redone on any cluster.
6. **The data bundle is about 22.1 GB and 3,296 files, 6% of `s2g0321`.** `s2g0321` holds 390 GB and 48,838 files in total. Every npz in the 9 dirs belongs to a kept pair, so filtering by pair saves nothing. Rorqual already has all of it.
7. **Two concurrency or inode bugs carry over from the census script:**
   - Depression-band seeds go to a shared `$SCR/seeds/<seed>_a00h12.json`. s5 and s6 of one model write the same file at the same moment. Write them per run instead.
   - Fit outputs and checkpoints are written into `rho_redesign/`, on rorqual's /project. All run outputs belong on /scratch.
8. **CUDA driver risk on other clusters.** The pip CUDA 12.9 NVVM produces PTX that needs a recent driver (CUDA ≥ 12.9 in `nvidia-smi`). Any cluster with an older driver fails at the first kernel launch. The smoke run on each cluster must print the driver version.

---

## A. How kimchi works (checked against the source)

### A.1 Where it runs for glusyn

| item | value (checked) |
|---|---|
| binary | `~/.local/bin/kimchi` → `/lustre09/project/6070394/dhuruva/kimchi/.venv/bin/kimchi` (same tree as `/project/rrg-emuller/dhuruva/kimchi`) |
| ledger | `/lustre09/project/6070394/dhuruva/kimchi_ledger_glusyn` (`kimchi.toml`, `events.jsonl`, `results/`, git repo, no `[ledger] remote`) |
| local cluster | rorqual (`local = name == $CC_CLUSTER`) |
| rorqual spool (glusyn) | `/lustre09/project/6070394/dhuruva/.kimchi_spool_glusyn` (on /project; move it to /scratch, see F) |
| other spool on rorqual | `/project/rrg-emuller/dhuruva/kimchi_spool` belongs to **nibi's** kimchi (MMCR/kettle); `claims.tsv` was updated today. Do not touch it. |
| plugin | `kimchi/kimchi/plugins/glusyn.py` → symlink to `glusynapse-fits/kimchi_plugin.py` (untracked in the kimchi repo, so a fresh kimchi clone needs it again) |
| ssh | `~/.ssh/config`: `Host nibi fir killarney vulcan tamia`, `HostName %h.alliancecan.ca`, `ControlPath ~/.ssh/cm/%r@%h:%p`, `ControlPersist yes`. The sockets work only on the login node that opened them. |

### A.2 The race and the claim

The rules are at the top of `sync.py`, in `decide_run()`:

- **One live copy everywhere.** An unclaimed, unfinished run keeps one live copy on every eligible cluster. A cluster is eligible when:
  - it offers one of the sweep's `gpu` keys;
  - it has a `dirs` entry for the project;
  - its observation succeeded;
  - it is on the agreed commit with a clean tree;
  - it is not excluded for that run;
  - fewer than 2 submit errors have been recorded for it.
- **The first copy to start wins.** The claim goes to the copy with the earliest `started` record (epoch seconds), from `<spool>/inbox`.
  - Other copies: pending ones are cancelled with `scancel --state=PENDING`. Running ones are cancelled only while younger than `grace_seconds` (600 s in our toml).
  - `act()` writes `claims.tsv` to every spool *before* it cancels or submits anything.
- **In-job check.** In the job, `job.sh single` looks the run up in `claims.tsv`. If another copy holds it, the job records `skipped` and exits at once. So a loser runs at most until the next sync, which comes every 60 s while anything is queued.
- **Failures.**
  - A copy that ended without a complete result (`is_complete` false) counts as an attempt. CANCELLED never counts. After `max_attempts` (2) the run is `failed`.
  - `OUT_OF_MEMORY` excludes that cluster for that run.
- **UNKNOWN clusters.** A cluster whose `observe` call fails is UNKNOWN: nothing on it is cancelled and nothing is submitted to it.
- **Who starts work.** `sync` and `status` never submit. Only `submit` does; `submit --watch` and `watch` loop. `watch` = `submit --watch` without the push.

### A.3 Code transport

`push_code()` in `sync.py` and `code()` in `agent.py`:

1. **Push.** On the local checkout: `git push <remote> HEAD:refs/heads/<branch>`.
2. **Fast-forward each other cluster.** On each cluster with a `dirs` entry: `git fetch <remote> <branch>`. If the checkout is on a different branch, it checks out that branch first. Then `git merge --ff-only FETCH_HEAD`. It never resets or stashes. A failure leaves that cluster out of the race.
3. **Commit and clean check.** `code_ok()` requires every reachable cluster for the sweep to be on one commit, or on `commit = "<sha>"` if the sweep pins one; `commit = "any"` turns the check off. A dirty tree is dropped unless `--allow-dirty` is given.
   - **Dirty means tracked files only**: `git status --porcelain -uno`. Untracked files, such as `results/` and a `results` symlink, do not make a tree dirty.
   - A *tracked* file modified at run time does make it dirty: `edge_params.npz` on a cache miss (see G.6).
4. **`--dirty`.** This pushes a snapshot commit of the working tree to `refs/kimchi/dirty`, and the other clusters check it out detached. Useful for tests. Not for the census, because the recorded commit would then be an unreachable snapshot.
5. **Remote clusters need GitHub read access.** Each remote cluster runs `git fetch` from GitHub itself, under BatchMode. Each cluster needs a GitHub-authorised key (a user key or a read-only deploy key).

### A.4 Project plugin format (`kimchi/plugins/__init__.py`, `sync.py`)

| hook | required | what kimchi does with it |
|---|---|---|
| `command(run_name, overrides, project_dir) -> (script, args, env)` | yes | Each manifest line is `key, workdir, env, script, args`. `job.sh` runs `cd workdir && env KIMCHI_RUN=<sweep>/<run> <env> bash <script> <args>`. |
| `results_path(project_dir, run_name, jobidraw) -> str` | yes | The agent reads this JSON on the winning cluster's **login node** (it must be on a shared filesystem). |
| `is_complete(results) -> bool` | yes | False means the copy counts as a failed attempt. |
| `KEEP = [...]` | no | **Top-level keys only.** The agent drops every other key before sending. Nested values are kept whole. |
| `metric(results) -> float` | no | The number shown in `status` and `results`. |
| `check(runs, project_dir) -> {run: error}` | no | Runs in-process on the login node before submit. Must be stdlib only and light. |
| `summarize`, `DERIVE` | no | `DERIVE[name](parent_results) -> {key: value}` for `after`/`derive` dependents. The parent and child may run on different clusters. |

The kimchi agent (`agent.py`) runs as `python3 -` on each cluster's login node over ssh. It is stdlib only and does git, squeue, sacct, file reads and sbatch. This is light and not a compute job. Nothing is installed remotely.

### A.5 Sweep files (`spec.py`)

- **Where they are read.** From `<local checkout>/sweeps/*.toml`, because `[projects.glusyn] sweeps = "sweeps"`. Only the top level is read; files starting with `_` are bases, and `retired/` is ignored.
- **Keys.** `name, project, gpu, time, run_seconds, name_template, race, hold, clusters, commit, max_attempts, exclude, common, axes, variants, extends, drop, retired, invalidate_on_commit, mem`.
- **The run grid.** The cartesian product of `[axes]` × `[variants.*]`. Each run's overrides are `[common]` + the variant's `overrides` + the axis values, rendered as `KEY=value` strings. Lists are joined with `,`; tabs and newlines are forbidden. They reach the script as positional args.
- **Run names.** `name_template.format(variant=<tag>, **axes)`. Names must match `[A-Za-z0-9][A-Za-z0-9_.+-]*` and be unique.
- **`extends`** merges `[common]`, `[axes]` and each variant key by key. `drop = ["common.X"]` removes an inherited key.
- **Dependents.** `after = "<variant>"` waits for the run with the same axes in that variant; `derive` adds overrides from the parent's results.
- **Archive.** Every submitted sweep is archived as `sweeps/archive/<name>.<digest>.json` in the ledger. Each result points to that archive.
- **Changing overrides.** Editing the overrides of a run that has already been submitted turns it `changed` and stops it. Rename the run instead.

### A.6 GPU keys, time buckets, mem, CPUs

- **GPU keys.** In `kimchi.toml`: `[clusters.X] gpus = { <key> = "<sbatch flags for ONE gpu>" }`. A sweep's `gpu` is a preference list; each cluster uses the first key it has (`gpu_for()`). A cluster with none of the keys is not raced.
- **Time buckets.** `{bucket}` in the flags becomes `b1` (≤3 h), `b2` (≤12 h), `b3` (≤24 h), `b4` (≤3 d) or `b5` from the run's time (`slurm.bucket`). Alliance routes partitions itself, so our flags do not need it.
- **Time.** `time` is a string, or a per-gpu dict. A cluster's `min_time` raises it (tamia 1 h).
- **Memory.** `mem` (sweep-level, or **per variant**: `[variants.X] mem = "82G"`) is appended as ` --mem=<mem>` **after** the cluster flags. sbatch takes the last value, so it overrides any `--mem` in the key. Runs with different mem go into separate arrays.
- **CPUs.** There is no sweep or variant key for CPUs. They come only from the cluster's GPU-key flags, so the glusyn keys must carry `--cpus-per-task=1`.
- **The sbatch call.**
  ```
  sbatch --parsable --account=<cluster account> <key flags> [--mem=..] --time=<t> \
    --job-name=k.<sweep>.<cluster>.<t>.<n> --output=<spool>/logs/%x_%A_%a.out \
    --export=ALL,KIMCHI_CLUSTER=..,KIMCHI_SPOOL=.. --array=0-<n-1> <spool>/job-<hash>.sh single <spool>/batches/<b>.tsv
  ```
  Arrays hold at most 500 tasks.
- **Job environment.** `--export=ALL` passes the environment of the `bash --login` shell the agent ran in, so `$SCRATCH` and `CC_CLUSTER` are set in the job.
- **Queue cap.** `max_queued = 800` per cluster counts all of the user's jobs (array tasks individually). The Alliance cap is 1000.
- **Packs.** `pack` (whole-node) clusters are opt-in. Per-run `mem` is multiplied by `gpus × oversubscribe` (`scale_mem`), not by `gpus`. With the template's `oversubscribe = 2`, a 92G run would ask for 736G on a 4-GPU node and fail. Use `oversubscribe = 1`, or leave tamia out (recommended, see D).

### A.7 Results

- **Inbox.** `job.sh` writes `started`, `finished rc` and `skipped` records to `<spool>/inbox`. `sync` reads them, then `fetch` reads `results_path` (KEEP-filtered) and stores `results/<sweep>/<run>.json` in the ledger.
- **What a result holds:**
  - from kimchi: cluster, copy, node, `start_t`/`end_t`, `commit`, `gpu`, `overrides`, `sweep_file`, `results_path`;
  - from the plugin: `"results": {KEEP keys}`.
- **Other copies.** Their results go to `results/<sweep>/attempts/`.
- **Logs.** `<spool>/logs/k.<batch>_<jobid>_<task>.out` on the cluster that ran the copy. They never come back, except a `log_tail` (last 30 lines) when the JSON is missing.

### A.8 Corrections to KIMCHI.md

| KIMCHI.md says | correct state |
|---|---|
| kimchi runs only on the nibi login node and is not installed on rorqual | Two installs. nibi runs MMCR/kettle (shared ledger). **rorqual runs glusyn** (own ledger `kimchi_ledger_glusyn`, rorqual local, smoke 22105657 done). |
| Draft repo `/lustre09/.../glusynapse-fits_draft/` | It is `/lustre09/project/6070394/dhuruva/glusynapse-fits`, with origin `git@github.com:dhuruvapriyan/glusynapse-fits.git`. Local `b6426e1` is not pushed; remote is at `8555fbb`. |
| Fit = `fit_v2.py --backend jax`, closure fit_v2/batch_v2/model_v2/jax_v2/targets | Stage 4 = `fit_launch.py` + `gpu_v10_rho`/`gpu_v11_rho` + `fit_v7` (numba-cuda, no jax). The closure is in C. |
| `h100_1g` (1g.10gb) slice, `--mem=45G`, `--time=03:30:00` (later 07:30) | Stage 4 uses **3g.40gb** (sjob.sh `-g`). Fits need 29G/63G/82G (stage-3 MaxRSS × 1.25) and val needs 92G. Fit time is about 25 min (extrapolated); see F.5. |
| Data = 5 `*_delta-prefire-vca` dirs (7.5 GB) + `basis_results_edges_sabrina_n120_delta` | Data = 9 split2 dirs under `/scratch/dhuruva/s2g0321/extracted` + `cexp` + 3 basis dirs + 3 keep lists (22.1 GB, 3,296 files), plus small csv/json in git. |
| Pins from `plastyfire/.venv` (jax 0.7.1, numba 0.67) | Stage 4 runs in `glusynapse_v2/.venv_v3` (`requirements_v3.txt`: numba 0.65.1, numba-cuda 0.30.2, cuda-toolkit 12.9, jax 0.7.1 also present). `env_v3.sh` unsets `CUDA_HOME`, so no cuda module is needed. |
| "CPUs cannot be set per sweep" | Correct. Also: `mem` can be set per variant, which KIMCHI.md does not mention. |
| Clean tree required | Tracked files only (`-uno`); `--dirty` snapshot mode exists. |
| `edge_params.npz` 47 KB covers every synapse | Stale. plastyfire's copy is 62 KB (updated 2026-10-02 13:52, with L2/3 synapses). Repo `targets.py` has also drifted from plastyfire. |
| `kimchi watch` must not be run bare (shared ledger) | True on nibi. On rorqual the glusyn ledger is separate, so a bare `watch` is safe there. |

---

## B. Data manifest

These are the inputs a stage-4 **fit** and its **val** read, from `run_stage4_fit.sh` lines 44–126 and `fit_v7.build` / `join_cexp` / `batch_v2`. `W=/scratch/dhuruva/s2g0321`, `T=delta-split2-ljp25g0321-prefire`. Both modes read the same data; the fit additionally reads the target csvs and the seed jsons.

### B.1 Large data (Globus or tar; not git)

| path | du -sh | files | read by |
|---|---|---|---|
| `$W/extracted/ebner_$T-vca` | 1.4G | 218 npz | `--dirs` (L5) |
| `$W/extracted/markram_$T-vca` | 943M | 152 | `--dirs` |
| `$W/extracted/sj03_$T-vca` | 2.6G | 197 | `--dirs` |
| `$W/extracted/sj03r50_$T-vca` | 750M | 44 | `--dirs` |
| `$W/extracted/sj07_$T-vca` | 792M | 85 | `--dirs` |
| `$W/extracted/l5extra_$T-vca` | 1.5G | 247 | `--dirs` (paired_l5_extra) |
| `$W/extracted/ebner_l23l5_$T-vseg-rs` | 5.6G | 438 | `--l23-dirs`, `--extra l23` |
| `$W/extracted/zilberter_l23l23_$T-vseg-rs` | 8.3G | 1,529 | `--extra l23l23` |
| `$W/extracted/l23l23extra_$T-vseg-rs` | 199M | 118 | `--extra l23l23` |
| `$W/cexp/{L5L5,L23L5,L23L23}.csv` | 525K | 3 (skip `cexp/parts/`) | `CEXP_DIR` → `join_cexp` (every basis synapse must be present) |
| `$W/basis_l5l5` | 649K | 22 csv | `ANALYTICAL_BASIS_DIR` |
| `$W/basis_l23l5` | 3.4M | 120 | `L23_BASIS_DIR` |
| `$W/basis_l23l23` | 3.4M | 120 | `--extra l23l23` basis |
| `$W/doublets/{l5,l23l5,l23l23}_keep_pairs.txt` | 3 KB | 3 (of 16 in the dir) | `K5`, `K235`, `K2323` |
| **total** | **≈ 22.1 GB** | **3,296** | |

The rest of `$W` is 368 GB and 45.5K files (cache, configs, logs, refitting_results, valid_*, `cexp/parts`, `cexp_sab`). None of it is read.

Npz names are `<pre>-<post>__<protocol>.npz`. I checked the counts against the keep lists: every npz in the 9 dirs belongs to a kept pair (22 L5 pairs = `subset24_pairs.txt` ∩ `l5_keep_pairs.txt`; 102 L2/3→L5; 120 L2/3→L2/3). Pair filtering would remove nothing.

### B.2 Small inputs (git)

| path (relative to plastyfire) | size | read by |
|---|---|---|
| `ebner/pair_geometry_L23PC_L5TTPC.csv` | 20 KB | `GEOM_L23` (L2/3 pair set + distal flags) |
| `ebner/ebner_targets.csv` | 15 KB | `targets.py` |
| `analytical_method/constants.py` | 2 KB | `batch_v2` |
| `nevian/nevian_protocols.py` | 20 KB | `targets.py` (lazy) |
| `glusynapse_v2/edge_params.npz` | 62 KB | `model_v2.edge_params` (the cache must cover every synapse; the h5 is not copied) |
| `glusynapse_v2/subset24_pairs.txt` | <1 KB | L5 pair set |
| `rho_redesign/r1D_s.csv`, `s1C_s_val_l23.csv`, `s1C_s_val_l23l23.csv` | 3 + 1 + 1 KB | fit mode: `dropof` builds the drop lists |
| seed jsons `rho_redesign/{s1C_s,s2E_s,s2C_s,s3W_s,s3C_u,s4V_u,s4V_s,s4N_u,s4Nm8_s}.json` + `/scratch/dhuruva/s1c_l23fit/s1CL_s.json` | 5–9 KB each | `SEEDJ` / `SEEDX` (the census uses the six s3*/s4* files) |

Not needed remotely:
- `data/dhuruva_modified_edges.h5` (4.3 GB): set `PLASTYFIRE_EDGES=/nonexistent`, so a cache miss fails loudly;
- `rho_redesign/results/v4_LM0chk.json` (only for `--check-v4`);
- the circuit and the simulation outputs.

### B.3 What to copy where, and the lean bundle

- **rorqual:** nothing. Point `GSF_DATA` at `/scratch/dhuruva/s2g0321`.
- **Each other cluster:** the 14 rows of B.1 (22.1 GB, 3,296 files), into **`$SCRATCH/gsf_data/s2g0321/`**, keeping the same relative layout. That way the run script only swaps the root (`W=$GSF_DATA`).
- **How to send it.** Two options:
  - **Recommended:** a Globus transfer of the 14 paths (directory sync) straight to `$SCRATCH` on the target. 3.3K files is nothing for Globus. The only inode-limited filesystem, rorqual's /project at about 467K/500K, is not touched, because the source is rorqual's /scratch.
  - **Tarball:** one uncompressed tar (npz is already compressed), `/scratch/dhuruva/gsf_bundle/s2g0321_stage4.tar` (≈22 GB, 1 inode). Keep it on each target's /project as a purge backup (1 inode there), and extract to /scratch. Build it in a small job, not on the login node:
    ```
    cd /scratch/dhuruva && tar -cf gsf_bundle/s2g0321_stage4.tar \
      s2g0321/extracted/{ebner,markram,sj03,sj03r50,sj07,l5extra}_delta-split2-ljp25g0321-prefire-vca \
      s2g0321/extracted/ebner_l23l5_delta-split2-ljp25g0321-prefire-vseg-rs \
      s2g0321/extracted/{zilberter_l23l23,l23l23extra}_delta-split2-ljp25g0321-prefire-vseg-rs \
      s2g0321/cexp/{L5L5,L23L5,L23L23}.csv s2g0321/basis_l5l5 s2g0321/basis_l23l5 s2g0321/basis_l23l23 \
      s2g0321/doublets/{l5,l23l5,l23l23}_keep_pairs.txt
    ```
    Sizing: no comparable job exists, so this is a pilot at 1 CPU, 1G, 0:30. Record `seff`.
- **Do not stage per job into `$SLURM_TMPDIR`.** 22 GB × 256 runs = 5.6 TB of reads.

---

## C. Code manifest

### C.1 Import and file closure of a stage-4 fit and val

From `fit_launch.py`, `gpu_v1{0,1}_rho.py`, `fit_v7.py` and recursively:

- **rho_redesign:** `fit_launch.py`, `fit_v7.py`, `fit_v6.py`. The kernel does `import fit_v6`; `fit_launch` aliases it to `$FITTER`, and `fit_v6.py` is needed when `FITTER` is unset.
- **Kernels:** `gpu_v10_rho.py`, `gpu_v11_rho.py`, `gpu_v7x_rho.py`, `gpu_v7_rho.py`, `gpu_v6_rho.py`, `gpu_v5_rho.py`, `gpu_v4_rho.py`, `rho_v4.py`.
- **glusynapse_v2:** `fit_v2.py`, `fit_v3.py` (`rules_row`), `gpu_v3.py` (imported by `gpu_v4_rho`), `batch_v2.py`, `model_v2.py`, `targets.py`, `edge_params.npz`, `subset24_pairs.txt`. `jax_v2.py` is imported only by `fit_v2 --backend jax`; keep it for the old sweeps.
- **Outside glusynapse_v2:** `analytical_method/constants.py`, `ebner/ebner_targets.csv`, `ebner/pair_geometry_L23PC_L5TTPC.csv`, `nevian/nevian_protocols.py`.
- **Drivers:** `rho_redesign/run_stage4_fit.sh` and `rho_redesign/stage4_env.sh`. `run_stage3_fit.sh` is only *referenced* in comments; nothing is sourced from it. Also `env_v3.sh`.
- **Shell tools in the job:** `jq`, `awk`, `comm`, `paste`, `nvidia-smi`. `jq` comes from the CVMFS StdEnv/2023 gentoo prefix, so it is present wherever CVMFS is.

The code totals 355 KB, about 25 files. Keep **plastyfire's relative layout** in the repo (`glusynapse_v2/…`, `glusynapse_v2/rho_redesign/…`, `ebner/…`). The `HERE`/`ROOT` path logic then works unchanged, and `diff -r` against plastyfire shows any drift.

### C.2 Repo state against what stage 4 needs

| in repo (`b6426e1`) | status |
|---|---|
| `glusynapse_v2/{batch_v2,fit_v2,jax_v2,model_v2}.py`, `subset24_pairs.txt`, `analytical_method/constants.py`, `ebner/ebner_targets.csv`, `nevian/nevian_protocols.py` | identical to plastyfire (`cmp`); keep |
| `glusynapse_v2/targets.py` | **stale** (differs from plastyfire); refresh |
| `glusynapse_v2/edge_params.npz` | **stale** (47 KB vs 62 KB); refresh |
| `requirements_v3.txt` | matches the `.venv_v3` pins; keep (it is the stage-4 env) |
| `env/bootstrap_cluster.sh` | useful; change it (venv v3 only, data on `$SCRATCH`, extra prints), see F |
| `env/setup_cluster.sh`, `requirements.txt` (jax env) | older path; keep for td4 sweeps or retire |
| `scripts/run_fit.sh`, `kimchi_plugin.py`, `sweeps/{_base,smoke,td4_pl5sj07}.toml`, `seeds/…td3_s5.json`, `fit_results/delta-cooker.json` | fit_v2/td4 era; move sweeps to `sweeps/retired/`. The plugin must be rewritten for stage 4 (the ledger's `glusyn` project uses this file). |
| `GLOBUS.md` | lists the old 7.5 GB `-prefire-vca` data; replace with B.1 |

**Missing:**
- `glusynapse_v2/fit_v3.py`, `glusynapse_v2/gpu_v3.py`;
- the whole `glusynapse_v2/rho_redesign/` closure from C.1;
- `ebner/pair_geometry_L23PC_L5TTPC.csv`;
- the three target csvs;
- the 10 seed jsons;
- a stage-4 run script, an env script, a summary writer, the new plugin and the new sweeps.

### C.3 Hard-coded paths that must become env-configurable

The default must stay the current rorqual value, so that plastyfire behaviour is unchanged. File:line refers to plastyfire as of today.

| file:line | hard-coded | proposed env (default = current value) |
|---|---|---|
| `rho_redesign/run_stage4_fit.sh:44` | `R=/project/rrg-emuller/dhuruva/plastyfire` | `GSF_ROOT` (repo root; `R=${GSF_ROOT:-…}`) |
| `rho_redesign/run_stage4_fit.sh:44` | `W=/scratch/dhuruva/s2g0321` (→ X, cexp, basis, doublets on lines 111, 113–117, 125) | `GSF_DATA` |
| `rho_redesign/run_stage4_fit.sh:45` | `SCR=/scratch/dhuruva/stage4` (phase A, seeds, lopo) | `GSF_SCR` (kimchi: `$GSF_OUT/<run>`) |
| `rho_redesign/run_stage4_fit.sh:54,55` | `/scratch/dhuruva/s1c_l23fit/s1CL_s.json` | `$GSF_ROOT/…/seeds/s1CL_s.json` (copy into git) |
| `rho_redesign/run_stage4_fit.sh:73` | `OUT=$RS/${NAME}_$S` (outputs inside the code tree, on /project) | `GSF_OUT` (+ run tag) |
| `rho_redesign/run_stage4_fit.sh:116` | `GEOM_L23` default `$R/ebner/…` | already env; follows `R` |
| `rho_redesign/run_stage4_fit.sh:207` | `o=$SCR/seeds/$(basename $j .json)_a00h…json` (shared name, concurrent writers) | write to the per-run dir |
| `env_v3.sh:5` | `source glusynapse_v2/.venv_v3/bin/activate` (relative to plastyfire) | `GSF_VENV` |
| `env_v3.sh:8` | `NUMBA_CACHE_DIR=/lustre09/project/6070394/dhuruva/plastyfire/glusynapse_v2/.numba_cache_v3` | `NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-$SCRATCH/.numba_cache_glusyn}` |
| `model_v2.py:38` | `EDGES` default `/project/rrg-emuller/.../dhuruva_modified_edges.h5` | already `PLASTYFIRE_EDGES`; set `/nonexistent` remotely |
| `rho_redesign/fit_v7.py:120`, `fit_v6.py:110` | `STANDIN_DIR = "/scratch/dhuruva/split1/cexp_standin"` | `GSF_STANDIN_DIR` (only reached with `CEXP_DIR=STANDIN`) |
| `rho_redesign/fit_v7.py:121`, `fit_v6.py:111` | `RELABEL_DIR = "/scratch/dhuruva/split1/cexp_relabel"`. **Written on every seeded start** (`relabel()` runs when `CEXP_DIR` is set and the first seed is unscaled). | `GSF_RELABEL_DIR` (kimchi: the per-run dir) |
| `rho_redesign/gpu_v7_rho.py:204`, `gpu_v7x_rho.py:167`, `gpu_v10_rho.py:732,764`, `gpu_v11_rho.py:927,959` | `CEXP_DIR` default `/scratch/dhuruva/split1/cexp` (`setdefault` / `get`) | harmless while the driver exports `CEXP_DIR`; leave |
| `sjob.sh:19,33,48` | `ROOT=/project/rrg-emuller/...`, account `def-emuller`, 3g.40gb gres, `.venv` | not used under kimchi: kimchi.toml supplies the flags |
| `rho_redesign/fit_v7.py:78,80`, `batch_v2.py:35` | ROOT-relative defaults for `GEOM_L23`, `L23_BASIS_DIR`, `ANALYTICAL_BASIS_DIR` | already env-overridable; the driver sets them |

Comment-only mentions (`fit_v7.py:18,28`, `fit_v6.py:8,18`, `gpu_v7_rho.py:14`, `run_stage4_fit.sh:28`) can stay. `submit_stage4_census.sh` and `stage4_env.sh` have no absolute paths. Under kimchi, the census script is replaced by a sweep file and `stage4_env.sh` by `[common]` overrides.

---

## D. Per-cluster table

"verify" means it cannot be checked from rorqual without ssh. `bootstrap_cluster.sh` should print it (F.3). The env is the same everywhere: `module load StdEnv/2023 python/3.11` + a venv from `requirements_v3.txt`. There is no cuda module, `CUDA_HOME` is unset, and the CUDA 12.9 runtime and NVVM come from pip wheels.

| cluster | GPU for us | MIG | account (kimchi.toml / docs) | checkout root (ledger comment / runbook) | scratch | RAM fit (29–82G fit, 92G val) |
|---|---|---|---|---|---|---|
| **rorqual** (local) | H100 80GB | **yes**: 1g.10gb/2g.20gb/3g.40gb (our jobs use `nvidia_h100_80gb_hbm3_3g.40gb`) | `def-emuller` (works in sjob.sh) | `/lustre09/project/6070394/dhuruva/glusynapse-fits` (exists) | `/scratch/dhuruva` (data already there) | **yes**: the census 82G/92G jobs on 3g slices are accepted |
| nibi | H100 80GB | verify (`sinfo -h -o %G`) | `def-emuller_gpu` (verify `sacctmgr`) | `/project/6055588/dhuruva/glusynapse-fits` | `/scratch/dhuruva` (verify) | verify; likely yes on 8-GPU nodes |
| fir | H100 80GB | verify | `def-emuller_gpu` (verify) | `/scratch/dhuruva/glusynapse-fits` (runbook; scratch purge risk for code + venv; prefer a /project path if one exists, verify) | `/scratch/dhuruva` | verify; likely yes |
| vulcan | L40S 48GB | **no** (L40S has no MIG) | `aip-emuller` | `/project/aip-emuller/dhuruva/glusynapse-fits` | verify | verify: 92G per GPU may exceed the per-GPU share of a 4×L40S node |
| killarney | L40S 48GB and H100 80GB | verify (H100 nodes) | `aip-emuller` | `/project/6102315/dhuruva/glusynapse-fits` | verify | verify |
| tamia | 4× H100 whole node (`pack`, opt-in, `min_time` 1 h) | no | `aip-emuller` | `/project/6102315/dhuruva` (same path as killarney in the runbook; verify the filesystem) | verify | **leave out**: whole nodes, the pack mem×2 issue (A.6), and a 25-min run on a full H100 is poor value |

The PAICE clusters (vulcan, killarney, tamia) bill `aip-emuller`. That this allocation may be used for this project is still the user's open decision (KIMCHI.md Q6). A whole L40S per fit costs more fairshare than a 3g slice (KIMCHI.md Q3, still open).

---

## E. The user's manual steps (in order)

An agent cannot do these: Duo, GitHub UI, Globus UI, and remote logins. Run everything on **rorqual3**, the login node that holds the live sessions, unless another cluster is named.

1. **Decide** (anywhere):
   - (a) PAICE (`aip-emuller`) allowed for glusyn: yes/no.
   - (b) Whole L40S lanes: yes/no.
   - (c) tamia out (recommended).
   - (d) The selection-rule details in H: HIT = |z| ≤ 1? The LOPO χ² definition? The tie tolerance? The held-out list?
2. **GitHub** (browser):
   - Confirm `github.com/dhuruvapriyan/glusynapse-fits` exists and is private. It exists: `git ls-remote` returns `8555fbb`.
   - For each new cluster, add that cluster's public key (step 4) as a **read-only deploy key** on the repo, or as a user SSH key.
3. **Open the sessions** (rorqual3, in tmux):
   ```
   kimchi login nibi fir killarney vulcan
   ```
   Approve one Duo push per cluster. fir, killarney and vulcan may answer `session live`. A login node reboot means doing this again.
4. **On each new cluster** (from rorqual3: `ssh nibi`, then `ssh fir`, …; the session from step 3 means no new Duo prompt):
   ```
   ls ~/.ssh/id_ed25519.pub || ssh-keygen -t ed25519 -N ""      # then add the .pub as in step 2
   ssh -T git@github.com                                         # expect "successfully authenticated"
   git clone git@github.com:dhuruvapriyan/glusynapse-fits.git <ROOT>/glusynapse-fits    # <ROOT> from table D
   cd <ROOT>/glusynapse-fits && bash env/bootstrap_cluster.sh <ROOT>                    # after the agent updates it (F.3)
   ```
   Then copy the printed block (gres list, `sacctmgr` accounts, `diskusage_report`, `which jq`, the paths in `~/.gsf_env`) back to the orchestrator.
   - The bootstrap builds the venv with `pip` on the login node. This is an install, not a computation.
5. **Globus** (web UI):
   - Source: the rorqual collection (verify the name, e.g. `alliancecan#rorqual`), `/scratch/dhuruva/s2g0321/…`, the 14 paths of B.1.
   - Destination: `$SCRATCH/gsf_data/s2g0321/` on each new cluster, same relative layout.
   - Option: transfer the single tar instead and extract it on the target in a small job:
     ```
     sbatch --account=<acct> --cpus-per-task=1 --mem=1G --time=0:30:00 --wrap "tar -xf <tar> -C $SCRATCH/gsf_data"
     ```
     This is a pilot size; check `seff`.
6. **Symlink results** (each new cluster, after step 4):
   ```
   mkdir -p $SCRATCH/glusyn_runs && ln -s $SCRATCH/glusyn_runs <ROOT>/glusynapse-fits/results
   ```
   This keeps outputs off /project. Untracked symlinks do not make the tree dirty.
7. **Smoke runs** (rorqual3, after F.5): `kimchi submit s4smoke_<cluster> --dry-run`, then `kimchi submit s4smoke_<cluster>`. Do one cluster at a time. Then run `kimchi submit s4smoke_<cluster> --watch` in tmux, or let the orchestrator drive it with a single `run_in_background` call.
8. **Census** (rorqual3, after the smoke runs pass and the heavy run is logged in DECISIONS.md): `kimchi submit s4census --dry-run`, then `kimchi submit s4census --watch` in tmux.

rorqual itself needs none of steps 2–6. The rorqual lane can start as soon as the agent has pushed (F.1–F.4).

---

## F. The implementing agent's steps

Rules: no python on the login node; jobs sized from measurements; stop after submitting and report job IDs. All repo edits go to `/lustre09/project/6070394/dhuruva/glusynapse-fits` (branch `main`). Do not edit plastyfire's driver; port it.

**F.0 Measure first.**
- When the running census 22314719–42 finishes, run `seff` on every `f4_*` and `v4_*` job.
- Record fit elapsed and MaxRSS per model (3W/4V/4N) and val per model in `run_stage4_fit.sh` and DECISIONS.md. Today these exist only as extrapolations (fit about 25 min).

**F.1 Code into the repo** (copy, keep the plastyfire layout; `scripts/sync_from_plastyfire.sh` gets the new file list):
- add `glusynapse_v2/fit_v3.py`, `glusynapse_v2/gpu_v3.py`;
- add `glusynapse_v2/rho_redesign/{fit_launch,fit_v7,fit_v6,gpu_v4_rho,gpu_v5_rho,gpu_v6_rho,gpu_v7_rho,gpu_v7x_rho,gpu_v10_rho,gpu_v11_rho,rho_v4}.py`;
- refresh `glusynapse_v2/targets.py` and `glusynapse_v2/edge_params.npz`;
- add `ebner/pair_geometry_L23PC_L5TTPC.csv`;
- add `glusynapse_v2/rho_redesign/{r1D_s.csv,s1C_s_val_l23.csv,s1C_s_val_l23l23.csv}`;
- add `seeds/{s1C_s,s2E_s,s2C_s,s3W_s,s3C_u,s4V_u,s4V_s,s4N_u,s4Nm8_s,s1CL_s}.json`;
- make the C.3 path edits in the **repo copies** (`fit_v7.py:121`/`fit_v6.py:111` `GSF_RELABEL_DIR`, `:120`/`:110` `GSF_STANDIN_DIR`), with defaults equal to the current values. Port the same edits back to plastyfire later, by hand.

**F.2 New scripts:**
- `env/env_v3.sh` sets:
  - `source ${GSF_VENV}/bin/activate`;
  - `unset CUDA_HOME`;
  - `NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-$SCRATCH/.numba_cache_glusyn}`;
  - `PLASTYFIRE_EDGES=${PLASTYFIRE_EDGES:-/nonexistent}`.
- `scripts/run_stage4.sh`, the kimchi entry point. It is a port of `run_stage4_fit.sh`:
  - `source ~/.gsf_env`, then export the `KEY=value` args;
  - `R=$GSF_ROOT`, `W=$GSF_DATA`, run dir `D=$GSF_ROOT/results/$GSF_RUN_TAG` (the results symlink to scratch);
  - `SCR=$D`, `OUT=$D/fit`, `SA=$D/fitA`, seeds in `$D/seeds`, `GSF_RELABEL_DIR=$D/relabel`;
  - generic rule knobs so that 64 variants need no script edits: `SETX` (json merged into `SET` after the MODEL case), `FREEX` (free params added), `KERNEL`, `SEEDX`, `TSET` (`keep5|k11|all` target set);
  - `VAL=1` runs the val mode on `$OUT.json` in the same allocation;
  - `HELDOUT=selection/heldout.txt` is appended to `--drop-targets` in both fit and val;
  - `VALJ=<json>` is a val-only rescore of a committed json, for the cross-GPU repro test;
  - log `nvidia-smi --query-gpu=name,driver_version --format=csv` and the `nvidia-smi` CUDA version;
  - at the end, read peak RSS from `/sys/fs/cgroup/memory.peak`, or from `sstat -j $SLURM_JOB_ID.batch -o MaxRSS` (verify which works per cluster);
  - write `results/<run>.json` with the summary writer.
- `scripts/summarize_stage4.py`, run **in the job** on the compute node. It builds the H.2 schema from `fit.json`, `fit*.csv`, `fit_val.json`, `fit_val*.csv` and `selection/heldout.txt`.
- `selection/heldout.txt` (the user's list, committed *before* the census) and `scripts/select_stage4.sh` (jq over the ledger results; H.3).
- `kimchi_plugin.py`:
  - `SCRIPT = "scripts/run_stage4.sh"`;
  - `KNOWN = {MODEL, START, KERNEL, SETX, FREEX, TSET, SEEDX, POP, MAXITER, MAXB, SIGMA, A00, MW, DEMOTE, FIXANC, HINGE, INITADM, STRATEGY, FITTER, VAL, VALJ, NAMEX, SMOKE}`;
  - `results_path = f"{project_dir}/results/{run_name}.json"`;
  - `KEEP` = the top-level keys of H.2;
  - `is_complete = de_fun is not None and chi2_total is not None` (or `val_status` for VALJ runs);
  - `metric = de_fun`.
- Update `env/bootstrap_cluster.sh`:
  - build only `gsf_venv_v3` from `requirements_v3.txt`;
  - write `~/.gsf_env` with `GSF_DATA=$SCRATCH/gsf_data/s2g0321`, `GSF_VENV`, `NUMBA_CACHE_DIR`;
  - `mkdir $SCRATCH/.kimchi_spool_glusyn`;
  - print `sinfo -h -o %G | tr , '\n' | sort -u`, `sacctmgr -nP show assoc user=$USER format=account`, `diskusage_report`, `which jq`.
- Rewrite `GLOBUS.md` from B.1. Update `.gitignore` (`results`, without the trailing slash, also covers the symlink). Move the old sweeps to `sweeps/retired/`.
- `diff -r` the closure against plastyfire, commit and push (the user allows pushes). `git status --porcelain -uno` must be empty.

**F.3 Ledger `kimchi.toml`** (`/lustre09/project/6070394/dhuruva/kimchi_ledger_glusyn/kimchi.toml`, glusyn only):
```toml
[clusters.rorqual]
account = "def-emuller"
spool = "/scratch/dhuruva/.kimchi_spool_glusyn"         # moved off /project (inodes); no live copies on the old one
gpus = { h100_3g = "--gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1 --cpus-per-task=1",
         h100_1g = "--gpus-per-node=nvidia_h100_80gb_hbm3_1g.10gb:1 --cpus-per-task=1" }   # 1g kept for the retired td4 sweeps

# per remote cluster, only once its bootstrap output shows the gres and the account:
# [clusters.nibi]  spool = "/scratch/dhuruva/.kimchi_spool_glusyn"
#   gpus = { h100_3g = "--gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1 --cpus-per-task=1" }   # only if that gres exists
# [clusters.vulcan] / [clusters.killarney] (only if the user accepted L40S)
#   gpus = { l40s = "--gpus-per-node=l40s:1 --cpus-per-task=1" }                               # gres name from sinfo

[projects.glusyn]
plugin = "kimchi.plugins.glusyn"
dirs = { rorqual = "/lustre09/project/6070394/dhuruva/glusynapse-fits" }   # + each bootstrapped cluster's checkout
remote = "origin"
branch = "main"
sweeps = "sweeps"
```
Do not add a full-H100 key unless the user asks for one.

**F.4 rorqual env:**
- `~/.gsf_env` on rorqual:
  - `GSF_DATA=/scratch/dhuruva/s2g0321`;
  - `GSF_VENV=/lustre09/project/6070394/dhuruva/plastyfire/glusynapse_v2/.venv_v3` (reuse it; saves 11.6K inodes);
  - `NUMBA_CACHE_DIR=/scratch/dhuruva/.numba_cache_glusyn`.
- `mkdir -p /scratch/dhuruva/glusyn_runs && ln -s /scratch/dhuruva/glusyn_runs /lustre09/project/6070394/dhuruva/glusynapse-fits/results`. The checkout's existing `results/smoke1` dir has to be moved into it first.
- A unit check of `summarize_stage4.py` on an existing stage-3 json goes in a small job (1 CPU, 2G, 0:15), never on the login node.

**F.5 Smoke sweeps**, one per cluster, `sweeps/s4smoke_<cluster>.toml`:
```toml
# 4N s5 at POP 15, 5+2 generations, then val, then the cross-GPU repro rescore.
# Sizing: smoke 22307548 18:14 / 57.7 GB; val 22300527 9:52 / 72.9 GB -> 28 min x 1.5 = 0:45, 72.9 x 1.25 = 92G.
name = "s4smoke_rorqual"
project = "glusyn"
gpu = ["h100_3g", "l40s"]
time = "00:45:00"
mem = "92G"
clusters = ["rorqual"]
name_template = "smk_{variant}_rorqual"
[common]
FITTER = "fit_v7"
INITADM = 5
HINGE = '{"l5/10Hz_10ms|control": [1.0, 100.0], "l5/10Hz_-10ms|control": [1.0, 100.0]}'
START = "s5"
POP = 15
MAXITER = 5
MAXB = 2
VAL = 1
[variants.fit4N]
tag = "4N"
overrides = { MODEL = "4N", SEEDX = "seeds/s4N_u.json,seeds/s4Nm8_s.json" }
[variants.repro]
tag = "repro"
overrides = { MODEL = "4N", VALJ = "seeds/s4N_u.json" }
```
Pass criteria:
- both runs `done`;
- the JSON shows `gpu_name`, a driver with CUDA ≥ 12.9 and `max_rss_gb`;
- `repro` `val_chi2.total` equals rorqual's to within 1e-6 on every cluster;
- the fit log shows no `skipped` pairs.

Submit the rorqual smoke run first. It needs no manual step after F.1–F.4.

**F.6 Census sweep**, `sweeps/s4census.toml` (generated from a variants table by a small bash/awk script when there are 64 variants):
```toml
name = "s4census"
project = "glusyn"
gpu = ["h100_3g", "l40s"]        # l40s only if accepted
time = "01:00:00"                # fit (F.0 seff, about 25 min est.) + val 9:52, x 1.5, rounded up; fix from F.0 before submitting
mem = "92G"                      # val peak 72.9 GB x 1.25; set per variant once the smoke JSON gives max_rss per model
name_template = "{variant}_{START}"
[common]
FITTER = "fit_v7"
INITADM = 5
HINGE = '{"l5/10Hz_10ms|control": [1.0, 100.0], "l5/10Hz_-10ms|control": [1.0, 100.0]}'
VAL = 1
[axes]
START = ["s5", "s6", "u7", "u8"]
[variants.m3W]
tag = "3W"
overrides = { MODEL = "3W", SEEDX = "seeds/s3W_s.json,seeds/s3C_u.json" }
[variants.m4V]
tag = "4V"
overrides = { MODEL = "4V", SEEDX = "seeds/s4V_u.json,seeds/s4V_s.json" }
[variants.m4N]
tag = "4N"
overrides = { MODEL = "4N", SEEDX = "seeds/s4N_u.json,seeds/s4Nm8_s.json" }
# ... rule variants as MODEL + SETX/FREEX/KERNEL/TSET overrides
```
Before submitting, log the heavy run in DECISIONS.md: 256 runs × 1 CPU × 1:00 = 256 CPU·h requested, about 150 3g-slice-h expected. Report the job IDs and stop; the orchestrator waits.

---

## G. Risks

1. **Fairshare per cluster.**
   - Each cluster's `def-emuller` / `aip-emuller` fairshare is separate, so spreading the load is the point.
   - Cancelled pending copies cost nothing. A loser that starts within about one poll (60 s while anything is queued) of the winner runs until the next sync, then gets cancelled (inside `grace_seconds` 600). The expected waste is 1–2 min per near-simultaneous start.
   - A whole L40S, or a full H100 if added, is billed well above a 3g slice for the same fit.
   - On MIG nodes the 92G host-RAM request may count against the slice's share of node RAM. Check `sshare` / the usage reports after the smoke run.
   - Run the smoke on each cluster before the census. Repeated sbatch errors (2 submit errors) exclude a cluster per sweep.
2. **Inodes.**
   - rorqual /project (about 467K/500K) gets no new files: data, venv and outputs are reused or kept on /scratch.
   - Move the spool to /scratch. It accumulates one log per started copy and one batch tsv per array.
   - Per remote cluster: the venv is about 11.6K files on that cluster's /project (its own quota; check `diskusage_report`); the data is 3.3K files on scratch; each run leaves about 21 files on scratch (fit, phase A, val, seeds, relabel, summary).
   - **Scratch purge** (60 days on the national clusters; verify on PAICE) can remove the data and the fir checkout. Keep the tar on /project as a backup.
3. **Outputs stay remote; only `results/<run>.json` returns.**
   - Everything the scorecard, the fig9 plots and the selection need must be inside that JSON (H.2), including the val rows and the full fit JSON.
   - The val **can and should** run in the same kimchi run (`VAL=1`). It then runs on the same GPU right after the fit, so its REPRO check is meaningful, and its numbers come back in the same JSON.
   - Cost: the run needs mem = max(fit, val), which is 92G even for 3W fits (29G alone). A separate `after` val variant could land on another cluster that lacks the fit JSON, so it is not worth it.
   - Large files (csv, ckpt) can be pulled later with Globus from `<scratch>/glusyn_runs/<run>/` on the cluster the ledger names.
   - `seff` on a remote job needs a login there, so the JSON carries `max_rss_gb` and elapsed for the sizing rule.
4. **Reproducibility across GPU types (L40S sm_89 against H100 sm_90).**
   - The kernels compute in float64 (float32 is only the stored inputs, upcast exactly) and use no atomics. So per-candidate objectives should agree to about 1e-12 across architectures; differences, if any, come from FMA contraction.
   - DE can still diverge after a near-tie selection. Starts on different GPUs are then different but valid samples. The basin check (≥ 2 starts within Δde_fun 2) covers this.
   - The `repro` smoke run measures it directly.
   - The pip CUDA 12.9 NVVM needs a driver that supports its PTX. Any cluster with an older driver fails the smoke run. Check the driver in the smoke JSON.
5. **Session fragility.** Masters live on rorqual3 only. A login-node reboot turns every remote cluster UNKNOWN: nothing is lost, but nothing is cancelled or submitted until `kimchi login` runs again. A running `watch` keeps working for the clusters still reachable.
6. **Silent-data traps:**
   - `batch_v2` skips pairs with no basis csv without an error. Check the "pairs:" and record counts in the smoke log, and add them to the JSON (`n_records` per model).
   - An `edge_params.npz` cache miss with `PLASTYFIRE_EDGES=/nonexistent` fails loudly, which is wanted. With a real h5 it would rewrite the tracked npz and make that checkout dirty.
   - A stale `seeds/*.json` in git versus plastyfire would silently change the seeding. Copy the seeds at F.1 and do not edit them in place.
7. **Concurrency bugs inherited from `run_stage4_fit.sh`.** The shared seed-copy names (C.3, line 207) and the shared `RELABEL_DIR` are fixed by the per-run dir in F.2.

---

## H. Selection-rule hook

### H.1 The approved rule, made mechanical

The rule is fixed before running: **Markram ±10 HIT first, then LOPO χ², then parameter count**, with a small held-out set that stays unseen until the end. These must be written down in DECISIONS.md **before** the census is submitted:

- **HIT:** both `l5/10Hz_10ms|control` and `l5/10Hz_-10ms|control` with |z| ≤ 1, z = (pred − mean)/SEM from the fit's unweighted binary-readout table. This matches the hinge's dead zone k = 1.0. Use `val_rows` to confirm.
- **Basin:** per variant, the best start by `de_fun`. The variant counts only if ≥ 2 of its 4 starts lie within Δde_fun ≤ 2 (as `pick` mode does today).
- **LOPO χ²:** run_stage4 fix 6 gives SE_pair per target from leave-one-L5-pair / leave-one-group rescores. The proposed definition is χ²_LOPO = Σ (pred − mean)² / (SEM² + SE_pair²) over the **core** targets. The target set and the definition need user sign-off.
- **Tie-break:** variants whose χ²_LOPO lies within a tolerance (e.g. 2) of the best are ordered by `n_free` (fewer first).
- **Held-out:** `selection/heldout.txt` is committed before the census. Those targets are in `--drop-targets` for both fit and val, so no run computes them. A final `heldout` rescore runs once, on the chosen model only.

LOPO is too costly to run on all 256 fits. Per fit: 22 L5 units × 2:31 + 10 L2/3→L5 groups × 5:00 + 10 L2/3→L2/3 groups × about 10 min, roughly 205 GPU-min (measured units 22307549/22307550 and val 22300527). Run it only for the best start of each HIT-passing variant, as a second sweep `s4lopo`:

- the agent writes those fit JSONs (from the ledger's `fit_json`) into `seeds/lopo/<run>.json`, commits and pushes;
- three runs per fit, one per pathway block, each aggregating its own SE_pair in-job and returning `lopo_chi2_<pathway>` + `sepair` rows.

Cheaper alternative (verify first): if each target prediction is a mean over pairs, the jackknife can be computed from per-pair predictions dumped by the val step, at no extra GPU cost.

### H.2 What every per-run JSON (`results/<run>.json`) must contain (= plugin KEEP)

| key | content | source | used for |
|---|---|---|---|
| `run`, `variant`, `model`, `start`, `seeded` | identity (`seeded` = START begins with s) | env | grouping, basin |
| `kernel`, `fitter`, `set`, `free`, `n_free` | rule definition; `n_free` = fitted dims (4 a's + free filters + gamma_d/p) | fit json `args`, `n_free` | parameter count (rule step 3) |
| `de_fun`, `de_fun_phaseA`, `nfev`, `minutes` | fitted weighted objective (binary readout) | fit json, `fitA.json` | basin pick, kimchi `metric` |
| `chi2`, `chi2_l23`, `chi2_l23l23`, `chi2_total`, `n_targets`, `n_targets_l23`, `n_targets_l23l23` | unweighted χ² per pathway on the fitted (core) targets | fit json | scorecard |
| `chi2_weighted`, `hinge_pen`, `chi2_ltd`, `chi2_ltp`, `aic`, `bic` | diagnostics | fit json | scorecard |
| `a`, `pre`, `gamma_d`, `gamma_p`, `v5`, `v10`, `v11` | fitted parameters | fit json | params table, reseeding |
| `markram` | `{"10Hz_10ms": {mean, sem, pred, z}, "10Hz_-10ms": {...}}` | fit `fit.csv` (l5) | **rule step 1** |
| `markram_hit` | bool, both |z| ≤ 1 | derived | rule step 1 |
| `val_chi2`, `val_n` | `{l5, l23, l23l23, total}`, all non-held-out targets, unweighted, rho_sigma 0 | val json | validation scorecard |
| `val_status` | `ok` / `repro_diff_accepted` / `failed` | val exit code | completeness |
| `val_rows` | `[[pathway, target, condition, mean, sem, pred, z, role]]`, role = core/validation; held-out rows absent | val csvs + heldout.txt | fig9-style plots, scorecard, LOPO χ² |
| `heldout_n`, `heldout_sha256` | count and hash of the held-out list in force | heldout.txt | proves the held-out targets stayed out |
| `n_records`, `n_pairs` | per model, from the log | run script | silent-skip guard (G.6) |
| `gpu_name`, `driver`, `cuda`, `host` | `nvidia-smi` | run script | cross-GPU repro |
| `max_rss_gb`, `elapsed_fit_s`, `elapsed_val_s` | cgroup peak / sstat; timestamps | run script | Slurm sizing rule (no remote seff) |
| `fit_json` | the whole fit json (≈ 8 KB) | fit json | rerun val/LOPO on any cluster |
| `lopo_chi2`, `lopo_chi2_<pathway>`, `sepair` | only in `s4lopo` runs, keyed by `parent_run` | lopo block | rule step 2 |

kimchi adds `cluster`, `node`, `commit`, `start_t`, `end_t`, `overrides` and `sweep_file`. With these fields, `scripts/select_stage4.sh` can apply H.1 to `kimchi_ledger_glusyn/results/s4census/*.json` and `results/s4lopo/*.json` with jq alone, on the login node, without fetching any remote file.

### H.3 The mechanical pick, as a script outline

1. `jq` over `results/s4census/*.json`: keep `.results.markram_hit == true`.
2. Group by `variant`. Within each group, sort by `de_fun`. Accept the group if ≥ 2 starts lie within 2 of the best, and keep the best start.
3. Write those runs' `fit_json` to `seeds/lopo/`, then commit and submit `s4lopo`.
4. Join on `parent_run`. Sort by `lopo_chi2`; within the tolerance, sort by `n_free`.
5. Report the top variant. Only then commit `VALJ=<winner> HELDOUT=` for the single held-out rescore.
