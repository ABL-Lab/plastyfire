# /project rrg-emuller inode audit (2026-10-01)

Proposal only: nothing has been moved, tarred or deleted. Each action needs user approval and runs as a small Slurm job via `glusynapse_v2/run_quota_action.sh`. The default mode is a dry run. Nothing in this plan runs on the login node.

**Quota.** rrg-emuller holds 496K of 500K files; space is fine (11 of 100 TB). /scratch/dhuruva has 4.2K of 1000K files and /nearline/rrg-emuller 9 of 5000 files.

**Source.** Scan job 22127461 (`run_quota_scan.sh`; 2:16, MaxRSS 818 MB). It wrote `/scratch/dhuruva/quota_scan/{dirs_depth1-5,hash,ext,caches,misc}.tsv` and `listing.tsv`. Our share is 333,957 inodes (files + dirs). That includes 4,159 ljp25 inodes, which are left out here because the offload job 22127439 is moving them.

**Never touched by this plan:**
- data/dhuruva_modified_edges.h5
- DEES_cell_packages
- `.git` internals (only `git gc` is proposed)
- inputs of active jobs:
  - 22126783 dendfit and 22127513/14 dfast: SSCxEModelExamples-abl-own/delta_caspike and plastyfire/.venv
  - 22127439 ljp25 offload
  - 22127440 nulls extraction: raw `delta-prefire-vseg-rs-nulls`
  - 22127441/42 fit_v4: `extracted/*_delta-prefire-vca` and `extracted/ebner_l23l5_delta-prefire-vseg-rs`

**In-use raw hashes, kept:**
- `delta-prefire-vseg` (Sabrina sims). The L5 `-vca` npz are extracted from it (`run_extract_vca_full.sh`).
- `delta-prefire-vseg-rs` and `delta-prefire-vseg-rs-nulls` (L2/3->L5).

## 1. Where the inodes are (scan, ljp25 excluded)

| subtree | inodes |
|---|---|
| plastyfire | 196,951 |
| - refitting_results/.../seed20262009 | 128,685 |
| -- Sabrina_L5TTPC_L5TTPC_STDP | ~113K (51,860 shared inputs + per-hash outputs below) |
| -- Ebner2019_L23PC_L5TTPC | ~18K (8,358 shared inputs) |
| - .venv (in use) | 32,674 |
| - glusynapse_v2 (.venv_v3 12,788 in use; extracted 9,301) | 23,700 |
| - analytical_method | 5,947 |
| - logs (843 in logs/nd_delta-cooker/) | 1,596 |
| MMCR_CL (.venv 55,910; .git 2,734; experiments/_rorqual_checkpoints 2,520) | 64,587 |
| abl-sabrina-paper (experiments/results/model_figures 19,598, git-ignored, active 2026-09-30) | 24,309 |
| plastyfitting (.venv 10,463; untouched since 2026-05-14) | 19,592 |
| SSCxEModelExamples-abl-own (input-traces 11,236 .ibw, git-tracked; used by running jobs) | 15,220 |
| DEES_cell_packages | 3,156 |
| kimchi (.venv 2,432) | 2,529 |

**Caches and build dirs** (`caches.tsv`):
- venvs: MMCR_CL 55.9K, plastyfire 45.5K (.venv + .venv_v3), plastyfitting 10.5K, kimchi 2.4K, top-level .venv 1.3K. Their pycache is included in these counts.
- `__pycache__` outside venvs: under 350 in total.
- `.ipynb_checkpoints`: 0. `node_modules`: 3. NEURON `x86_64`: 283 in plastyfire.
- `*.SUCCESS`: 840, all from the delta-cooker run.

**Top extensions** (`ext.tsv`):
- plastyfire: .h5 21K, .pkl 20K, .log 19K, .npy 18K, .json 17K, .npz 16K.
- abl-sabrina-paper: .pickle 12K.
- SSCx: .ibw 11K.

**Raw prefire outputs per hash** (`hash.tsv`). Each count includes the `bluecellulab_results_<hash>/` contents, `pool_<hash>_*.log` and `simulation_edges_<hash>.pkl`.

| hash | dataset | inodes | bytes | npz | status |
|---|---|---|---|---|---|
| delta-prefire | Sabrina 26,818 + Ebner 3,145 | 29,963 | 6.8 TB | partly (ebner_, sj03_, sj03r50_, ebner_l23l5_delta-prefire) | old baseline, superseded by vseg / vseg-rs |
| delta-prefire-tr | Sabrina | 5,872 | 134 GB | markram_delta-prefire-tr | rejected |
| antic-delta-prefire | Sabrina | 5,716 | 216 GB | analytical_method/extracted_sabrina_n120_antic-delta_d0p25 | rejected (antic-delta is still the emodel comparison model) |
| delta-sv-prefire | Sabrina | 4,704 | 1.33 TB | ebner_delta-sv, markram_delta-sv | rejected |
| optimizer | Sabrina | 4,200 | 217 GB | none | rejected/old |
| delta-sv-prefire-tr | Sabrina | 1,176 | 29 GB | none | rejected |
| delta-cooker neurodamus run | Sabrina + Ebner + logs | 15,659 | ~15 GB | compare_nd_bcl_delta-cooker.csv | rejected |
| delta-prefire-vseg | Sabrina | 5,390 | 1.18 TB | -vseg, -vca | **keep** |
| delta-prefire-vseg-rs / -nulls | Ebner | 3,145 / 1,590 | 0.35 TB | yes | **keep** |

The delta-cooker run's 15,659 inodes break down as:
- full_delta-cooker bcl, pool and pkl: 3,349
- pool_delta-cooker and simulation_edges_delta-cooker.pkl: 1,672
- neurodamus outputs (out/soma/rho.h5, *.dat, cooker json + .SUCCESS, pydamus logs): 8,400
- Ebner full_delta-cooker-full: 1,395
- plastyfire/logs/nd_delta-cooker: 843

**Empty dirs.** There are 5,239 empty `out/` dirs in refitting_results.

## 2. Git status (untracked bulk)

- **plastyfire:**
  - 58 untracked entries (scripts and pngs, plus dirs glusynapse_v2/, fit_results/, ebner/, basis_results_edges_* and threshold_results_*).
  - refitting_results and .venv are ignored.
  - glusynapse_v2 (23.7K) is untracked, so it has no git backup.
- **MMCR_CL:**
  - Untracked: NeurReps26-Overleaf-NSEE/ and toy_nsee/.
  - .venv is ignored; it is uv-managed and uv.lock is present.
  - 1,916 loose objects.
- **abl-sabrina-paper:**
  - 114 untracked entries: biodata/, experiments/data/branches/*.json and dees_cpost_values/.
  - experiments/results (20K) is ignored.
- **plastyfitting:**
  - Untracked: bluecellulab_results_fitted/ (4.3K) and basis_results_edges_mini/.
  - .venv is ignored; 242 loose objects.

## 3. Candidates

**How to submit.** `cd /project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2` and then:

```
sbatch -J q_<ITEM> --export=ALL,ITEM=<ITEM>,MODE=<MODE> [--time=H:MM:00] run_quota_action.sh
```

- Defaults: 1 CPU, 1G, 0:15, `--account=def-emuller`. Output goes to `/scratch/dhuruva/quota_scan/actions/`.
- **Always run `MODE=dryrun` first.** It writes the item list and its inode count to /scratch and changes nothing.
- **Tars:** `MODE=tar` writes one `$D/archive_tars/<ITEM>.tar`, with paths relative to $D (restore: `tar -xf <tar> -C /project/rrg-emuller/dhuruva`). `MODE=delete` refuses unless the tar entry count equals the inode count.

**Sizing.** No tar job has been measured yet. Memory is estimated at well under 1G, from scan 22127461 (818 MB, most of it awk tables). Time is estimated as 2 × bytes / 300 MB/s, covering the write plus the verify read. **Run `cooker` first as the pilot**, record its `seff` in the script, and rescale the other times from it.

**Tar location.** Tars stay on /project, one inode each. Optionally copy them to /nearline/rrg-emuller later (tars only). /scratch is not a good home for them: files are purged after ~60 days without access.

### Category (a): safe to delete

| item | path | inodes freed | why | job |
|---|---|---|---|---|
| A1 `empty_out` | `refitting_results/.../simulations/*/*/out/` (empty) | 5,239 | empty neurodamus output dirs; `rmdir` refuses non-empty dirs | `MODE=dryrun`, then `MODE=delete` (0:15) |
| A2 `sv` | delta-sv-prefire raw (Sabrina) | 4,704 | rejected; npz in extracted/ebner_delta-sv and markram_delta-sv | dryrun, then delete (0:15) |
| A3 `tr` | delta-prefire-tr raw (Sabrina) | 5,872 | rejected; npz in extracted/markram_delta-prefire-tr | dryrun, then delete (0:15) |
| A4 `gc` | MMCR_CL/.git, plastyfitting/.git | ~2,100 (est.) | packs 1,916 + 242 loose objects; no data lost | `ITEM=gc MODE=gc` (0:15) |

### Category (b): archive (tar, then delete originals after the counts match)

| item | path | inodes freed | bytes | --time | why |
|---|---|---|---|---|---|
| C1 `cooker` (pilot) | delta-cooker neurodamus run: per-pair outputs, pools and pkls (Sabrina), full_delta-cooker-full (Ebner), plastyfire/logs/nd_delta-cooker | 15,658 | ~15 GB | 0:15 | rejected, but the neurodamus results are costly to rerun |
| C2 `optimizer` | optimizer raw (Sabrina) | 4,199 | 217 GB | 0:45 | rejected/old; no npz |
| C3 `antic` | antic-delta-prefire raw (Sabrina) | 5,715 | 216 GB | 0:45 | rejected, but antic-delta is the emodel comparison model, so tar rather than delete |
| C4 `svtr` | delta-sv-prefire-tr raw (Sabrina) | 1,175 | 29 GB | 0:15 | rejected; no npz |
| C5 `plastyfitting` | all of plastyfitting except .venv (tarred, .git included); .venv deleted after `pip freeze` | 19,591 | 368 GB | 1:00 | untouched since 2026-05-14, superseded by plastyfire. Rebuild the venv from `/scratch/dhuruva/quota_scan/actions/plastyfitting_requirements.txt` |
| C6 `am_extracted` | plastyfire/analytical_method/extracted* | 5,910 | 33 GB | 0:15 | old analytical-method npz, superseded by glusynapse_v2/extracted. One dir has group `dhuruva` |
| C7 `v2x_rejected` | glusynapse_v2/extracted/{ebner_delta-sv, markram_delta-sv, markram_delta-cooker, markram_delta-prefire-tr, ebner_preview, *_delta-ljp25-prefire-vca} | 3,487 | 32 GB | 0:15 | npz of rejected variants (ljp25 rejected 2026-09-30). **Run only after the ljp25 offload 22127439 has finished.** |
| C8 `delta_prefire` | delta-prefire raw (Sabrina + Ebner) | 29,962 | 6.8 TB | ~13:00 if tarred | **needs the user's status decision.** It is not on the rejected list, but it is superseded and partly extracted. If confirmed obsolete, deleting it is much cheaper than a 6.8 TB tar |

### Relocate

| item | path | inodes freed | why | job |
|---|---|---|---|---|
| R1 `mmcr_venv` | MMCR_CL/.venv (7.1 GB) | 55,910 | a re-creatable venv. /home has room: 261K of 500K inodes, 27 of 50 GB. Copying avoids needing network access on a compute node | step 1: `ITEM=mmcr_venv MODE=copy --time=0:30:00 --mem=1G` copies it to `$HOME/venvs/mmcr_cl` and prints src/dst inode counts. Step 2, after approval and when no MMCR job is running: a 0:15 job with `--wrap "rm -rf /project/rrg-emuller/dhuruva/MMCR_CL/.venv && ln -s $HOME/venvs/mmcr_cl /project/rrg-emuller/dhuruva/MMCR_CL/.venv"`. The symlink keeps shebangs and `uv run` working. Check with the MMCR_CL session first; the project is active |

### Category (c): keep

- **In use:**
  - plastyfire/.venv (32.7K): running dendfit/dfast jobs and glusynapse_v2 scripts.
  - glusynapse_v2/.venv_v3 (12.8K): the default v3 fit kernel.
  - raw hashes delta-prefire-vseg, -vseg-rs and -vseg-rs-nulls (10.1K).
  - in-use extracted dirs.
- **Shared prefire inputs** (Sabrina 51.9K, Ebner 8.4K): needed to rerun any hash. They could be tarred per dataset once the fits are final.
- **abl-sabrina-paper/experiments/results/model_figures** (19.6K): active. Tar delta-cooker (9.8K) and antic-delta-cooker (8.9K) once the figures are frozen.
- **SSCxEModelExamples-abl-own** (15.2K): running jobs; input-traces are git-tracked.
- **DEES_cell_packages, data/dhuruva_modified_edges.h5:** never touched.
- **Small venvs:** kimchi/.venv 2.4K and top-level .venv 1.3K.

## 4. Ranked plan (lowest risk first)

| rank | item | inodes freed | cumulative | risk |
|---|---|---|---|---|
| 1 | A4 `gc` | ~2,100 | 2.1K | none |
| 2 | A1 `empty_out` | 5,239 | 7.3K | none (rmdir on empty dirs) |
| 3 | C1 `cooker` (pilot tar) | 15,658 | 23.0K | low; data kept in a tar |
| 4 | C4 `svtr`, C2 `optimizer`, C3 `antic` (tars) | 11,089 | 34.1K | low; data kept |
| 5 | C6 `am_extracted`, C7 `v2x_rejected` (tars; C7 after 22127439) | 9,397 | 43.5K | low; data kept |
| 6 | C5 `plastyfitting` (tar + venv drop) | 19,591 | 63.1K | low; the venv must be rebuilt if the project is ever reused |
| 7 | A2 `sv`, A3 `tr` (deletes, npz exist) | 10,576 | 73.7K | small; raw traces are lost, but they can be regenerated from the shared inputs |
| 8 | R1 `mmcr_venv` | 55,910 | **129.6K** | medium-low; brief MMCR downtime, must be coordinated |
| alt 8 | C8 `delta_prefire` (delete or tar) | 29,962 | 103.6K | needs the user's decision on the hash; a tar takes ~13 h |

Each tar adds one file, which is already subtracted above.

**Recommendation:** ranks 1-8 free about **130K inodes** (496K down to about 366K). Only ranks 2 and 7 delete data without keeping a tar, and in both cases the npz files or empty dirs mean nothing needed is lost. The ljp25 offload already running frees another ~4.2K.

**Order of operations:**
- Before each step, run its dryrun and check the inode count against this report.
- After each tar job, run `seff <id>` and record it in `run_quota_action.sh`.
- Re-run `run_quota_scan.sh` (1G, 0:15) at the end to confirm.
