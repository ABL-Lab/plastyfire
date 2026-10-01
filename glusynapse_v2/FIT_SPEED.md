# GPU DE fit speed (v2 → v3), 2026-09-30

pl5sj07 setup: 5 -vca dirs, paired_l5 + sjostrom07 (29 targets), 502 records, 3997 synapses, popsize 8 × 14 = 112.

## Why v2 was slow
1. **The time scan is bound by memory traffic.** `jax_v2` runs a `lax.scan` over up to 376k steps. Its carry is 8 float64 (P × lanes) arrays, and it is written to GPU memory and read back on every step. That is about 1.1e9 lane-steps × P × 128 B, or about 15 TB per call at P = 112. The cost is linear in P: 5.7 s at P = 5 and 79 s at P = 112 on 2g (22107148).
2. **The call grows because more candidates pass the rules.** Rules 1-4 send only admissible candidates to the GPU. That is 9/112 in the initial population and 112/112 at every checkpoint (nit 93-115). This took the 1g call from 13 s to 108-147 s. It is not a leak or cache growth.
3. **Recompiles.** Each new admissible count P is a new jit shape, so all 4 chunk shapes recompile. That costs about 6 s at P = 5 and about 10 s at P = 112 (22107148, rep 0 vs rep 1). It matters only in the early generations, because P = 112 is fixed later.
4. **Minor costs:**
   - 941 of the 4938 lanes are mglu/NO copies that re-run the whole rho/T/N/Z scan.
   - Bucket padding adds 1.12× work.
   - With CHUNK the traces are streamed to the GPU on every call (about 14 GB).
   - td/tp and the rules are Python loops over records, about 1.2 s per call.

## What v3 changes (new files; v2 untouched and kept as the reference)
- **`gpu_v3.py`** is one numba-cuda kernel with one thread per (candidate, synapse).
  - The whole time loop runs with its state in registers.
  - One thread carries the 3 dpre variants (control, A_mglu = 0, A_NO = 0), so there are no copy lanes.
  - Each synapse runs to its own length (no padding).
  - Arrivals and impulses are event lists, and filter coefficients are recomputed only when h changes.
  - Traces stay on the GPU (7.14 GB, fits a 1g.10gb slice). It compiles once for any P.
  - It uses the same float64 operations in the same order as v2. dt and the numerics are unchanged.
- **`fit_v3.py`** uses the same DE, rules, seeding and checkpoints as `fit_v2`, and `--resume` works across v2 and v3. td/tp for the whole population are computed once, and host traces are freed after the upload.
- **Environment:** `env_v3.sh` activates `glusynapse_v2/.venv_v3`, with its own numba cache directory.
  - The venv pins the same numpy 1.26.4, scipy, pandas and jax 0.7.1 as `plastyfire/.venv`.
  - It adds numba 0.65.1 + numba-cuda 0.30.2 from the Compute Canada wheelhouse. The built-in numba.cuda cannot target sm_90 with the pip NVVM.

## Equality (the test any rule change must pass in both implementations)
Tests: `tests/test_jax_v3.py` via `run_test_jax_v3.sh` (TDRIVE=4 or 3), and `tests/smoke_gpu_v3.py` for all t_drive/no_drive/pre_drive modes on synthetic data.

| test | job | χ² rel. diff | pred diff | rho / dpre diff | rho binary flips |
|---|---|---|---|---|---|
| td4, 5 fixed + 112 population | 22107148 | 1.2e-15 | 8.9e-16 | 1.8e-15 / 4.7e-15 | 0 |
| td3, 5 fixed + 112 population | 22107522 | 4.7e-16 | 4.4e-16 | 1.4e-15 / 7.3e-15 | 0 |
| smoke, 7 modes | 22107147 | – | – | 4.4e-16 max | – |

The two implementations agree to rounding, so the optimum cannot move.

## Speed-up (late generations, all candidates admissible)
| slice | v2 objective per generation | v3 per generation | speed-up |
|---|---|---|---|
| 1g.10gb | 108-147 s (22098444-7) | 4.5-4.9 s at P 91-95 (22107523) | ~25-30× |
| 2g.20gb | 79-80 s at P 112 (22107148) | 2.2-2.5 s at P 91-95 (22107524); 2.95 s at P 112 | ~27× |

- **Full fit:** v2 projected about 5 h for 150 generations on 1g (07:30 requested). The production resumes 22107773-6 on 1g ran their last 33-38 generations at 4.2-5.8 s/gen, in 3:57-4:46 wall including about 2 min to load the data, and reached χ² 39.45-41.79.
- **150-generation fit on v3:** about 2 min load + 150 × 5 s ≈ 15 min on 1g, or about 9 min on 2g.
- **Why not 100×:** the kernel is now bound by float64 arithmetic (about 9e10 synapse-steps per generation) on 1/7 of an H100. Going further would need float32 state, which changes the numerics near the rho threshold, or fewer steps. Neither was done.

## Sizing (seff)
| job | MaxRSS | wall | CPU eff. |
|---|---|---|---|
| resumes 22107773-6 | 21.6-26.3 GB | 3:57-4:46 | 89-98% |
| pilots 22107523/4 | 26.3 / 23.8 GB | – | – |
| equality tests (v2 + v3 models in one process) | 40.1 GB | 8:13 | – |

For a full v3 fit: 1 CPU, `--mem=33G`, `--time=00:30:00` on 1g.10gb.

## How to run
See the header of `run_reduced_gpu_v3.sh`. It takes the same env vars as `run_reduced_gpu.sh`, and CHUNK is ignored.
