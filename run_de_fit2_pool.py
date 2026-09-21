"""
Run the 700 DE-FIT-2 (1 ms model) BCL STDP simulations as ONE Slurm job
with a process pool,
instead of 700 one-per-job sbatch submissions.

Why
---
Each simulation takes ~1m10s (max ~17m for pairs that miss the c_pre/c_post
cache and have to run a live threshold search), but submit_ic_best_700.sh asks
for --time=01:30:00 per job. With 700 jobs the queue wait, not the compute,
sets the wall clock: jobs submitted at 01:28 started at 03:24 to do 90 seconds
of work. Packing them into one 60-way pool turns ~700 x (2h queue + 1.2 min)
into a single queue wait plus ~700 * 1.2 / 60 ~= 15 minutes of compute.

Each task runs pairrunner_edges_fit.py as a subprocess with cwd=workdir,
exactly as the sbatch --chdir did, so the per-simulation behaviour is
unchanged — same script, same args, same output filename.

Usage (inside an salloc/sbatch allocation):
    python run_de_fit2_pool.py --workers 60
    python run_de_fit2_pool.py --workers 60 --dry-run
    python run_de_fit2_pool.py --workers 60 --skip-existing
    python run_de_fit2_pool.py --workers 60 --skip-existing --force
        # retry only the workdirs with no pkl, ignoring the induction
        # spike-count guardrail (see --force below)
"""

import argparse
import multiprocessing
import os
import subprocess
import sys
import time

PLASTYFIRE_ROOT = "/lustre06/project/6077694/dhuruva/plastyfire"
DEFAULT_STDP_NAME = "L5TTPC_L5TTPC_STDP"

def results_dir_for(stdp_name):
    return os.path.join(
        PLASTYFIRE_ROOT,
        f"refitting_results/fitting/n100/seed19091997/{stdp_name}/simulations",
    )

RESULTS_DIR = results_dir_for(DEFAULT_STDP_NAME)
PAIRRUNNER = os.path.join(PLASTYFIRE_ROOT, "plastyfire/pairrunner_edges_fit.py")
CACHE = os.path.join(PLASTYFIRE_ROOT, "cpre_cpost_cache/ion_channels_tau278.pkl")
CIRCUIT_CONFIG = os.path.join(
    PLASTYFIRE_ROOT, "data/dhuruva_modified_ion_channels_circuit_config.json"
)
# NOTE: non-canonical label. md5 of the *discarded* gamma_p=738.83 set;
# FIT_ARGS below were refit under the physiological gamma bounds without
# changing this string. The canonical md5 of the params actually used is
# 1f1909bc65dc. Kept as-is because 658 result pkls already carry this name.
PARAM_HASH = "b8c7ff3ecf0a"

PROTOCOLS = ["10Hz_10ms", "10Hz_-10ms", "10Hz_5ms", "10Hz_30ms",
             "10Hz_-30ms", "10Hz_50ms", "10Hz_-50ms"]

# DE fit #2: refit at 1 ms (stride 4 on extracted_d0p25) after dt_scan showed
# the 2 ms model carried maxcurve 0.01045 with 6 knife-edge synapses.
#
# Measured in bluecellulab with these params (hash b8c7ff3ecf0a, 658 pairs):
#
#     dt   bluecellulab   in vitro
#    -10     0.8065        0.7922
#     +5     1.2376        1.2038
#    +10     1.1549        1.2013
#
# (The analytical prediction and weighted err for this particular run were not
# written to a log; the numbers above are the ones actually measured.)
#
# Refit under the physiological gamma bounds gamma_d in [50, 200] and
# gamma_p in [150, 300]: the first attempt railed gamma_p at 738.83, outside
# range, and was discarded. The kept fit sits inside both bounds
# (gamma_d 77.76, gamma_p 299.91 -- the latter is at the ceiling, so treat any
# further potentiation demand as a sign the model is under-potentiating).
#
# hash b8c7ff3ecf0a. No threshold re-injection needed for bluecellulab
# (a-params are passed, so simulator_edges computes theta at runtime from
# cpre_cpost_cache). Neurodamus DOES need injection -- see
# submit_neurodamus_defit.py.
FIT_ARGS = [
    "--gamma_d_GB_GluSynapse=77.7558",
    "--gamma_p_GB_GluSynapse=299.9121",
    "--a00=1.003498", "--a01=2.902478",
    "--a10=1.644558", "--a11=2.764812",
    "--a20=1.003498", "--a21=2.902478",
    "--a30=1.644558", "--a31=2.764812",
    "--tau_effca_GB_GluSynapse=278.3177658387",
    f"--param_hash={PARAM_HASH}",
    f"--cpre-cpost-cache={CACHE}",
    f"--circuit-config={CIRCUIT_CONFIG}",
]

# ---------------------------------------------------------------------------
# c_pre-only method (--cpre-only)
# ---------------------------------------------------------------------------
# Reduced Graupner-Brunel thresholds: theta_d = a00*c_pre, theta_p = a10*c_pre,
# with the c_post coefficients pinned to zero (a01 = a11 = a21 = a31 = 0).
# Four free params instead of six: a00, a10, gamma_d, gamma_p.
#
# Fitted with  analytical_method/fit.py --cpre-only --fit-gamma --de
#              --anchor-tails --sign --workers 60
#
#   a00 = 1.004709   a10 = 1.977768   (apical tied: a20=a00, a30=a10)
#   gamma_d = 84.3845   gamma_p = 186.5408
#   weighted err 7.718 in 2592 evaluations  (vs 2.594 for DE fit #2)
#
# Offline analytical prediction for this set:
#
#     dt   analytical   in vitro
#    -50     0.8594        -
#    -30     0.8402        -
#    -10     0.8384      0.7922
#     +5     1.1612      1.2038
#    +10     1.1242      1.2013
#    +30     0.9885        -
#    +50     0.8883        -
#
# a01/a11/a21/a31 MUST be passed explicitly as 0.0, not omitted:
# simulator_edges._apply_theta_from_a_params does fit_params.get("a01", 1.0),
# so a missing key silently becomes a coefficient of 1.0.
#
# CPRE_ONLY_HASH is canonical -- it is the md5 of exactly these params, unlike
# DE fit #2's b8c7ff3ecf0a (see the note above).
CPRE_ONLY_HASH = "bdbf06f915d0"
CPRE_ONLY_ARGS = [
    "--gamma_d_GB_GluSynapse=84.3845",
    "--gamma_p_GB_GluSynapse=186.5408",
    "--a00=1.004709", "--a01=0.0",
    "--a10=1.977768", "--a11=0.0",
    "--a20=1.004709", "--a21=0.0",
    "--a30=1.977768", "--a31=0.0",
    "--tau_effca_GB_GluSynapse=278.3177658387",
    f"--param_hash={CPRE_ONLY_HASH}",
    f"--cpre-cpost-cache={CACHE}",
    f"--circuit-config={CIRCUIT_CONFIG}",
]


# --cpre-only-pm10: c_pre-only refit scoring ONLY dt=-10 and +10 (fit.py
# --fit-dt=-10,10, no tail anchors). User-supplied 2026-09-03:
#   a00 1.014763  a10 1.809655  gamma_d 161.8697  gamma_p 281.6517
#   weighted err 1.736 over 2 points / 3392 evals
#   offline curve: -50 0.8493  -30 0.8186  -10 0.8033  +5 1.1235
#                  +10 1.1233  +30 0.9880  +50 0.9080
# Only -10/+10 entered the objective, so +5 and the tails are UNSCORED
# predictions here, not fitted values.
# a01/a11/a21/a31 must be passed explicitly as 0.0: simulator_edges uses
# fit_params.get("a01", 1.0), so a MISSING key silently becomes coefficient 1.0.
CPRE_ONLY_PM10_HASH = "21f3bec952eb"   # canonical md5 of sorted(fit_params)
CPRE_ONLY_PM10_ARGS = [
    "--gamma_d_GB_GluSynapse=161.8697",
    "--gamma_p_GB_GluSynapse=281.6517",
    "--a00=1.014763", "--a01=0.0",
    "--a10=1.809655", "--a11=0.0",
    "--a20=1.014763", "--a21=0.0",
    "--a30=1.809655", "--a31=0.0",
    "--tau_effca_GB_GluSynapse=278.3177658387",
    f"--param_hash={CPRE_ONLY_PM10_HASH}",
    f"--cpre-cpost-cache={CACHE}",
    f"--circuit-config={CIRCUIT_CONFIG}",
]


# --cpre-only-s1: stride-1 (full 0.25 ms grid, no decimation) refit, the
# accurate-integration successor to CPRE_ONLY_PM10_ARGS. User-supplied
# 2026-09-07 from `fit.py --cpre-only --fit-gamma --de --seed 3 --popsize 24
# --tol 0.001 --stride 1 --fit-dt=-10,5,10 --anchor-tails`:
#   a00 1.020209  a10 1.810090  gamma_d 168.6454  gamma_p 298.1457
#   weighted err 5.673 / 19104 evals
# WARNING from fit.py itself: gamma_p is railed at 99% of [150,300] -- "the
# objective wants to leave the physiological range, so the fit is
# compensating for something else." NOT yet established as trustworthy;
# running through BCL specifically to check whether that railing produces a
# worse measured result the way the pm10 refit did (offline 1.74 -> BCL
# 10.15, worse than the original 4.70 -> 8.20). Compare against bdbf06f915d0
# (not railed, BCL 8.20/10.86) before treating this as an improvement.
# offline curve: -50 0.8809  -30 0.8422  -10 0.8194  +5 1.1650
#                +10 1.1103  +30 0.9902  +50 0.9108
CPRE_ONLY_S1_HASH = "283d79c04294"   # canonical md5 of sorted(fit_params)
CPRE_ONLY_S1_ARGS = [
    "--gamma_d_GB_GluSynapse=168.6454",
    "--gamma_p_GB_GluSynapse=298.1457",
    "--a00=1.020209", "--a01=0.0",
    "--a10=1.810090", "--a11=0.0",
    "--a20=1.020209", "--a21=0.0",
    "--a30=1.810090", "--a31=0.0",
    "--tau_effca_GB_GluSynapse=278.3177658387",
    f"--param_hash={CPRE_ONLY_S1_HASH}",
    f"--cpre-cpost-cache={CACHE}",
    f"--circuit-config={CIRCUIT_CONFIG}",
]


def find_tasks(skip_existing, results_dir=None, limit=None):
    """Return [(pair, protocol, workdir)] for every runnable workdir."""
    results_dir = results_dir or RESULTS_DIR
    tasks, skipped = [], 0
    for pair in sorted(os.listdir(results_dir)):
        pair_dir = os.path.join(results_dir, pair)
        if not os.path.isdir(pair_dir):
            continue
        for proto in PROTOCOLS:
            wd = os.path.join(pair_dir, proto)
            if not os.path.isfile(os.path.join(wd, "prefire_simulation_config.json")):
                skipped += 1
                continue
            if skip_existing and os.path.isfile(
                    os.path.join(wd, f"simulation_edges_{PARAM_HASH}.pkl")):
                skipped += 1
                continue
            tasks.append((pair, proto, wd))
            if limit and len(tasks) >= limit:
                return tasks, skipped
    return tasks, skipped


# Appended to every pairrunner invocation. Set from --force in main() BEFORE the
# pool forks, so the children inherit it (same copy-on-write trick the module
# globals rely on elsewhere in this file).
EXTRA_ARGS = []


def _worker(task):
    """Run one simulation in its own process, cwd=workdir (as sbatch --chdir did)."""
    pair, proto, wd = task
    t0 = time.time()
    log_path = os.path.join(wd, f"pool_{PARAM_HASH}_{proto}_{pair}.log")
    cmd = [sys.executable, PAIRRUNNER] + FIT_ARGS + EXTRA_ARGS
    try:
        with open(log_path, "w") as logf:
            subprocess.run(cmd, cwd=wd, stdout=logf, stderr=subprocess.STDOUT,
                           check=True)
        ok = os.path.isfile(os.path.join(wd, f"simulation_edges_{PARAM_HASH}.pkl"))
        return (pair, proto, ok, time.time() - t0,
                None if ok else "ran but no pkl written")
    except subprocess.CalledProcessError as e:
        return pair, proto, False, time.time() - t0, f"exit {e.returncode} (see {log_path})"
    except Exception as e:                                   # pragma: no cover
        return pair, proto, False, time.time() - t0, str(e)


def main():
    global FIT_ARGS, PARAM_HASH, EXTRA_ARGS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=60,
                    help="Parallel simulations (default: 60). Each peaks at ~4 GB, "
                         "so keep workers * 4 GB under the job's --mem.")
    ap.add_argument("--skip-existing", action="store_true",
                    help=f"Skip workdirs that already have simulation_edges_{PARAM_HASH}.pkl")
    ap.add_argument("--force", action="store_true",
                    help="Pass --force to pairrunner_edges_fit.py: ignore the "
                         "induction spike-count guardrail so workdirs that fired "
                         "51-60 post spikes instead of 50 still write a pkl "
                         "(marked guardrail_forced=True).")
    ap.add_argument("--cpre-only-pm10", action="store_true",
                    help="Run the c_pre-only set refit on dt=-10/+10 only "
                         f"(PARAM_HASH {CPRE_ONLY_PM10_HASH}).")
    ap.add_argument("--cpre-only-s1", action="store_true",
                    help="Run the stride-1 c_pre-only refit "
                         f"(PARAM_HASH {CPRE_ONLY_S1_HASH}). gamma_p is railed "
                         "at 99% of its bound -- see comment above CPRE_ONLY_S1_ARGS.")
    ap.add_argument("--cpre-only", action="store_true",
                    help="Run the reduced c_pre-only parameter set (theta_d = "
                         "a00*c_pre, theta_p = a10*c_pre, c_post coefficients "
                         "zero) instead of DE fit #2. Carries its own "
                         f"PARAM_HASH ({CPRE_ONLY_HASH}) so results never mix.")
    ap.add_argument("--stdp-name", default=DEFAULT_STDP_NAME,
                    help=f"*_STDP subtree under refitting_results/fitting/n100/"
                         f"seed19091997/ to run (default: {DEFAULT_STDP_NAME}). "
                         "e.g. L23PC_L5TTPC_STDP. NOTE: cpre_cpost_cache/"
                         "ion_channels_tau278.pkl only has entries for the "
                         "L5TTPC_L5TTPC_STDP pairs -- any other set is a 100% "
                         "cache miss, so simulator_edges falls back to live "
                         "c_pre/c_post mini-sims per task (2 extra short NEURON "
                         "sims per workdir, not cached across the 7 protocols "
                         "of a pair). That path is real but its time/MaxRSS are "
                         "unmeasured here -- use --limit for a pilot batch first.")
    ap.add_argument("--limit", type=int, default=None,
                    help="Only run the first N tasks found (for sizing a pilot "
                         "batch before committing --mem/--time on a full run).")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # Swap the parameter set BEFORE find_tasks() (which keys --skip-existing on
    # PARAM_HASH) and before the pool forks, so workers inherit it by
    # copy-on-write exactly as they do EXTRA_ARGS.
    n_methods = sum([args.cpre_only, args.cpre_only_pm10, args.cpre_only_s1])
    if n_methods > 1:
        sys.exit("--cpre-only, --cpre-only-pm10, --cpre-only-s1 are mutually exclusive")
    if args.cpre_only:
        FIT_ARGS, PARAM_HASH = CPRE_ONLY_ARGS, CPRE_ONLY_HASH
    elif args.cpre_only_pm10:
        FIT_ARGS, PARAM_HASH = CPRE_ONLY_PM10_ARGS, CPRE_ONLY_PM10_HASH
    elif args.cpre_only_s1:
        FIT_ARGS, PARAM_HASH = CPRE_ONLY_S1_ARGS, CPRE_ONLY_S1_HASH

    if args.force:
        EXTRA_ARGS = ["--force"]

    results_dir = results_dir_for(args.stdp_name)
    if not os.path.isdir(results_dir):
        sys.exit(f"no such results dir: {results_dir}")
    tasks, skipped = find_tasks(args.skip_existing, results_dir, args.limit)
    print(f"stdp set : {args.stdp_name}")
    print(f"tasks    : {len(tasks)}  (skipped {skipped})")
    print(f"workers  : {args.workers}")
    print(f"cache    : {CACHE}")
    print(f"method   : {'cpre-only stride-1 (gamma_p railed)' if args.cpre_only_s1 else 'cpre-only pm10 (fit on -10/+10 only)' if args.cpre_only_pm10 else 'cpre-only (a01=a11=0)' if args.cpre_only else 'DE fit #2'}")
    print(f"hash     : {PARAM_HASH}")
    print(f"force    : {args.force}")
    if args.dry_run:
        for pair, proto, _ in tasks[:10]:
            print(f"  [dry-run] {pair} / {proto}")
        print(f"  … {len(tasks)} total")
        return

    if not tasks:
        print("Nothing to do.")
        return

    t_start = time.time()
    done = failed = 0
    failures = []
    # maxtasksperchild=1: each simulation gets a fresh interpreter so NEURON/HOC
    # state cannot leak between pairs (the same reason the serial path forks a
    # subprocess per pair).
    with multiprocessing.Pool(args.workers, maxtasksperchild=1) as pool:
        for pair, proto, ok, secs, err in pool.imap_unordered(_worker, tasks):
            if ok:
                done += 1
            else:
                failed += 1
                failures.append((pair, proto, err))
            n = done + failed
            rate = n / max(1e-9, time.time() - t_start)
            eta = (len(tasks) - n) / rate if rate > 0 else 0
            print(f"[{n}/{len(tasks)}] {pair}/{proto} "
                  f"{'ok' if ok else 'FAIL: ' + str(err)} ({secs:.0f}s)  "
                  f"eta {eta/60:.1f} min", flush=True)

    print("=" * 70)
    print(f"Done: {done} ok, {failed} failed in "
          f"{time.strftime('%H:%M:%S', time.gmtime(time.time() - t_start))}")
    if failures:
        print(f"Failures ({len(failures)}):")
        for pair, proto, err in failures[:40]:
            print(f"  {pair}/{proto}: {err}")
        if len(failures) > 40:
            print(f"  … and {len(failures) - 40} more")


if __name__ == "__main__":
    main()
