# Fit-procedure fixes: code spec (user-approved 2026-10-02, FIT_META_DIAG.md section 6)

Status: part 1 and part 2 written (fit_v7.py, fit_launch.py, 2026-10-02); tests 22307548 (fit_v6 reference) vs the fit_v7 REPRO job, and the all-flags smoke.
Files: `scorecard.sh` (fix 1), `run_stage3_fit.sh` FIXANC block (fix 3), `run_stage4_fit.sh` (fixes 2, 3, 5, 6 driver
side, plus DEMOTE=1 and the a00 > 1 seeds), this file (fixes 2 and 4 fitter side).

## Where things live (read-only references; line numbers as of 2026-10-02 19:00)

- The kernels (`gpu_v10_rho.py`, `gpu_v11_rho.py`) only patch and call `fit_v6.main()` (gpu_v11 l.951-1001). The DE,
  the objective, the admissibility rules and the seed handling are all in **`fit_v6.py`**:
  - objective `objective_gpu` l.202-224: `F3.rules_row(td, tp, peak)` (fit_v3.py l.24-35: 1e4 if theta_d > theta_p at
    any L5 synapse, 1e3 if < 30 % / 15 % of L5 synapse-records cross theta_d / theta_p) on the L5 model `_M[0]` only;
    admissible members go to the GPU in one batch.
  - weights `install_weights` l.238-250 (instance wrapper of the class chi2, only when `--weights` is given).
  - seeds l.448-476: `pre = {**seed_set, **fj["pre"]}`, so **`--seed-set` only fills keys missing from the json's
    `pre`; it cannot override them, and a00..a11 come from the json's `"a"` block**. `--seed-set '{"a00": 1.2}'` does
    nothing. Seeded runs: uniform random population (rng `--seed`), seeds in members 1..n, `--x0` a's in member 0.
    Unseeded: `init="latinhypercube"`. `--resume NPZ` replaces the population with the saved one.
  - DE call l.492-494: `differential_evolution(..., popsize=args.popsize (default 8), init=init, polish=False,
    seed=args.seed, vectorized=True, updating="deferred", tol=1e-6)`, strategy default best1bin.
- Readout: `rho_sigma` (rho_v4 default 0) is read per candidate in `gpu_v4_rho.GPUModelV4.chi2` l.324-338 and
  `_ratios` l.296-322 (ndtr((rho - 0.5) / sigma) instead of the step; frozen lanes keep the step). No subclass up to v11
  overrides `chi2` or `_ratios`, and `fit_v2.unpack` starts P from the filters, so `"rho_sigma": s` in `--set` switches
  the smooth readout on for the DE **and** for the csv tables. No kernel change is needed.

## Fix 2: multi-start DE from admissible starts

Driver (done, `run_stage4_fit.sh`): one job per START (`s5 s6` seeded, `u7 u8` unseeded; more with `u9`...), popsize 15,
600 generations; `run_stage4_fit.sh MODEL - pick` ranks the starts by `de_fun` (the fitted objective, never validation)
and accepts the basin when at least 2 starts are within df <= 2 (`RS/<NAME>_pick.txt`). Seeded starts also get the
depression-band copies of each seed (a00 -> 1.2, a10 shifted by the same amount, so theta_p - theta_d is unchanged at
every synapse; made with jq because `--seed-set` cannot do it).

Fitter (part 2), in a **new file `fit_v7.py` = copy of `fit_v6.py`** (fit_v6 stays untouched for the running jobs):
1. argparse, after l.328:
   `--init-admissible` int, default 0 (off: the code path is unchanged, REPRO bit-identical) = max sampling rounds;
   `--init-pool` int, default 20 (candidates per round = pool x NP); `--strategy` str, default `"best1bin"`.
2. new helper after `objective_gpu` (l.224):
   ```python
   def rules_pen(X, chunk=512):
       """rules_row penalty for each row of X (n, D), on the L5 model, in GPU chunks."""
       out = np.empty(len(X))
       for s in range(0, len(X), chunk):
           AP = [unpack5(x) for x in X[s:s + chunk]]
           td, tp = _M[0][0].thetas_all([a for a, _ in AP], [{**_CFG["filters"], **P} for _, P in AP])
           out[s:s + len(AP)] = [F3.rules_row(td[i], tp[i], _CFG["peak"]) for i in range(len(AP))]
       return out
   ```
3. replace l.468-476 (init block) by: build `pop` as today (seeded: uniform + seeds in 1..n; unseeded: a Latin
   hypercube `lo + qmc.LatinHypercube(d=len(bounds), seed=rng).random(NP) * (hi - lo)` with
   `rng = np.random.default_rng(args.seed)`, `NP = args.popsize * len(bounds)`). Then, if `args.init_admissible > 0`:
   keep the seed rows; for every other row with `rules_pen > 0`, replace it by the next admissible candidate from
   LHS rounds of `pool x NP` points (at most `init_admissible` rounds); if still short, keep the lowest-penalty
   candidates and print a warning. Log one line:
   `admissible init: <a> of <NP> admissible (seeds <s>, rounds <r>, acceptance <a/cand>)`. Always set `init = pop`.
   (The unseeded default stays `init="latinhypercube"` when the flag is off.)
4. DE call l.492: add `strategy=args.strategy`. Diag suggestion: `rand1bin` for the first ~100 generations; with the
   two-phase driver that is `STRATEGY=rand1bin` in phase A only (one more env in run_stage4_fit.sh if wanted).
5. Driver: `FITTER=fit_v7 INITADM=5` (the driver already passes `--init-admissible` / `--strategy` when the env is set).
Check: generation-1 log "admissible" count must equal NP (today 3-16 of 56-72).

## Fix 4: Markram +-10 as a hard constraint (hinge)

Fitter (part 2, `fit_v7.py`):
1. argparse: `--hinge` default `os.environ.get("HINGE", "")`, json `{"l5/10Hz_10ms|control": [k, lam], ...}`;
   empty = off (REPRO). Global `HINGE = {}` next to `WEIGHTS` (l.227).
2. `install_weights` l.238-250: build `H = {k: HINGE[key] for k in T for key in (f"{k[0]}|{k[1]}", f"{name}/{k[0]}|{k[1]}") if key in HINGE}`
   and in `chi2w` replace the loop body by
   ```python
   z = (v - T[k][0]) / T[k][1]
   if k in H:                                   # hinge replaces z^2 and its weight for this target
       kk, lam = H[k]
       out += lam * np.where(np.isfinite(z), np.maximum(np.abs(z) - kk, 0.0) ** 2, 1e3)
   else:
       out += W[k] * np.where(np.isfinite(z), z ** 2, 0.0)
   ```
   (a NaN prediction of a hinge target costs 1e3 instead of 0, see FIT_META_DIAG 1c).
3. l.397: install the wrapper when `args.weights or args.hinge`; log `hinge <name>: {...}`.
4. json: `out["hinge"] = HINGE` and `out["hinge_pen"]` = the hinge part of the objective at the final point (recompute
   from the csv tables: sum lam max(|z| - k, 0)^2 over the hinge keys).
5. Values: `HINGE='{"l5/10Hz_10ms|control": [1.0, 100.0], "l5/10Hz_-10ms|control": [1.0, 100.0]}'`: zero inside +-1 SEM;
   1.5 SEM costs 25, 2 SEM costs 100 (the size of a whole stage-3 objective), so a miss is never traded away. The
   Markram weight MW is then irrelevant (the hinge ignores it). The csv chi2 stays unweighted z^2. The scorecard's
   Markram line checks the result.

## Fix 5: smoothed readout

Driver (done): phase A runs the DE with `"rho_sigma": SIGMA` (default 0.1) in `--set`; phase B `--resume`s phase A's
final population at rho_sigma 0 for MAXB (150) generations, so the reported json, csvs, de_fun and the val rescore are the
binary readout. SIGMA=0 runs one phase. No fitter change needed. An in-run anneal (change sigma inside the DE) is not
recommended: with `updating="deferred"` the stored population energies would be stale after each change; the restart in
phase B re-evaluates them.
Open: the sigma > 0 chi2 is a Python loop over targets (gpu_v4_rho l.329-337); the smoke job 22307548 measures the
generation time. If it is much slower, use a smaller MAXA or vectorise that loop in part 2 (new class in a new file).

## Fix 6: pair-sampling error from a leave-one-out rescore

Driver (done): `run_stage4_fit.sh MODEL START lopo IDX` (or an array): IDX 0-21 drop one L5 pair (L5-only rescore,
no extras), 22-31 drop one of 10 round-robin groups of the 102 L2/3->L5 pairs, 32-41 the same for the 120 L2/3->L2/3
pairs; `--maxiter 0` at the json in LOPOJ (REPRO DIFF accepted). `lopoagg` (pure awk) writes `RS/<TAG>_sepair.csv`
with SE_pair = sqrt((n-1)/n sum (p_i - mean)^2) per target and factor = SEM^2 / (SEM^2 + SE_pair^2).
`SEPAIR=<csv>` multiplies each target's DE weight by that factor, which is exactly chi2 with
SEM_eff^2 = SEM^2 + SE_pair^2 in the objective; the csv z stays on the data SEM, and the scorecard is unchanged.
SE_pair depends on the parameters a little: compute it at the current best json (s4N_u), fit, and recompute once at the
new best if the parameters moved a lot. No fitter change is needed. (Optional part-2 alternative: `--sem-extra CSV` in
`build` setting `T[k] = (m, sqrt(s^2 + se^2))`; not recommended, it changes the reported z.)

## Part 2: how to run fit_v7 without editing any gpu_v*.py

New file `RS/fit_launch.py` (the driver already calls it when env FITTER is set):
```python
"""Run a gpu_v1x kernel with fit_v6 replaced by the module in env FITTER (e.g. fit_v7)."""
import importlib, os, runpy, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.dirname(HERE)]
kernel = sys.argv.pop(1)
sys.modules["fit_v6"] = importlib.import_module(os.environ.get("FITTER", "fit_v7"))
sys.argv[0] = kernel
runpy.run_path(kernel, run_name="__main__")
```
The kernels' `import fit_v6` then returns fit_v7, and all their patches (`F.pack`, `unpack5`, `gpu_v6_rho`,
`BatchV2`, `MV`, `_CFG`, `_M`) land on it. Tests, in a Slurm job, 1 CPU, MIG slice: (a) `FITTER=fit_v7` with no new flag,
`--maxiter 0` on s4N_u.json and s3V_s.json: REPRO OK (bit-identical to fit_v6); (b) `INITADM=5`, MAXITER 2: log line
"admissible init: NP of NP"; (c) HINGE on, MAXITER 2: hinge line in the log, `hinge_pen` in the json.

## Not in this round
- Pathway-uniform admissibility rules (diag fix 8) and taking unreachable rows out of the core list (fix 5 of the diag)
  were not in the approved six.
- No new model parameters were added (the user allowed a few with anchors). The anchored values remove 1-3 free
  parameters (k 6 for every MODEL now).
