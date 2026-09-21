"""
Run the 700 DE-FITTED (analytical_method) BCL STDP simulations as ONE Slurm job
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
    python run_de_fit_pool.py --workers 60
    python run_de_fit_pool.py --workers 60 --dry-run
    python run_de_fit_pool.py --workers 60 --skip-existing
"""

import argparse
import multiprocessing
import os
import subprocess
import sys
import time

PLASTYFIRE_ROOT = "/lustre06/project/6077694/dhuruva/plastyfire"
RESULTS_DIR = os.path.join(
    PLASTYFIRE_ROOT,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations",
)
PAIRRUNNER = os.path.join(PLASTYFIRE_ROOT, "plastyfire/pairrunner_edges_fit.py")
CACHE = os.path.join(PLASTYFIRE_ROOT, "cpre_cpost_cache/ion_channels_tau278.pkl")
CIRCUIT_CONFIG = os.path.join(
    PLASTYFIRE_ROOT, "data/dhuruva_modified_ion_channels_circuit_config.json"
)
PARAM_HASH = "9e43ab966326"

PROTOCOLS = ["10Hz_10ms", "10Hz_-10ms", "10Hz_5ms", "10Hz_30ms",
             "10Hz_-30ms", "10Hz_50ms", "10Hz_-50ms"]

# DE fit against Markram 1997 (analytical_method/fit.py, 3888 evaluations,
# weighted err 2.594). Apical tied to basal. gamma_d/gamma_p were free too.
#
# Predicted by the offline model, which reproduces bluecellulab to
# max |EPSP ratio error| 0.0104:
#     dt   predicted   in vitro
#    -10      0.8112     0.7922
#     +5      1.1953     1.2038
#    +10      1.1791     1.2013
# Tails (not measured, anchored softly at 1.0) came out -50:0.889 -30:0.866
# +30:1.002 +50:0.888 -- still depressed, which is the main thing this BCL run
# is testing.
#
# hash 9e43ab966326. gamma_d/gamma_p are HOC globals and part of the hash, so
# this cannot collide with the tied-ic4 run (cdf3a1e1db98).
#
# No re-injection of thresholds is needed: a-params are passed, so
# simulator_edges takes the use_a_params branch and computes theta at runtime
# from cpre_cpost_cache. edges.h5 theta_d/theta_p are ignored on this path.
FIT_ARGS = [
    "--gamma_d_GB_GluSynapse=139.6809",
    "--gamma_p_GB_GluSynapse=483.4465",
    "--a00=1.100479", "--a01=2.891376",
    "--a10=1.543869", "--a11=3.095199",
    "--a20=1.100479", "--a21=2.891376",
    "--a30=1.543869", "--a31=3.095199",
    "--tau_effca_GB_GluSynapse=278.3177658387",
    f"--param_hash={PARAM_HASH}",
    f"--cpre-cpost-cache={CACHE}",
    f"--circuit-config={CIRCUIT_CONFIG}",
]


def find_tasks(skip_existing):
    """Return [(pair, protocol, workdir)] for every runnable workdir."""
    tasks, skipped = [], 0
    for pair in sorted(os.listdir(RESULTS_DIR)):
        pair_dir = os.path.join(RESULTS_DIR, pair)
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
    return tasks, skipped


def _worker(task):
    """Run one simulation in its own process, cwd=workdir (as sbatch --chdir did)."""
    pair, proto, wd = task
    t0 = time.time()
    log_path = os.path.join(wd, f"defit_pool_{proto}_{pair}.log")
    cmd = [sys.executable, PAIRRUNNER] + FIT_ARGS
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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=60,
                    help="Parallel simulations (default: 60). Each peaks at ~4 GB, "
                         "so keep workers * 4 GB under the job's --mem.")
    ap.add_argument("--skip-existing", action="store_true",
                    help=f"Skip workdirs that already have simulation_edges_{PARAM_HASH}.pkl")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    tasks, skipped = find_tasks(args.skip_existing)
    print(f"tasks    : {len(tasks)}  (skipped {skipped})")
    print(f"workers  : {args.workers}")
    print(f"cache    : {CACHE}")
    print(f"hash     : {PARAM_HASH}")
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
