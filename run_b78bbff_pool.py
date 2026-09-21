"""
Step 3 for the b78bbff arm: STDP induction at +/-10 ms, pooled into one job.

Mirrors run_ic_best_pool.py exactly -- same pairrunner, same gen-4 fit params --
but points --circuit-config at the b78bbff emodel and namespaces the output as
simulation_edges_b78bbff_<hash>.pkl so nothing collides with the production or
nseg results.

Unlike the nseg sweep this DOES pass --cpre-cpost-cache: here the emodel is the
variable, and its cache (cpre_cpost_cache/b78bbff_tau278.pkl) was computed on this
very emodel, so reading it is correct and avoids re-running the mini-sims 200 times.
"""

import argparse
import multiprocessing
import os
import re
import subprocess
import sys
import time

PLASTYFIRE = "/lustre06/project/6077694/dhuruva/plastyfire"
RESULTS_DIR = os.path.join(
    PLASTYFIRE,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations")
PAIRRUNNER = os.path.join(PLASTYFIRE, "plastyfire/pairrunner_edges_fit.py")
CACHE = os.path.join(PLASTYFIRE, "cpre_cpost_cache/b78bbff_tau278.pkl")
CIRCUIT = os.path.join(PLASTYFIRE, "data/dhuruva_b78bbff_circuit_config.json")
TAG = "b78bbff"
PARAM_HASH = "0ce64fa83b85"
PROTOCOLS = ["10Hz_10ms", "10Hz_-10ms"]

FIT_PARAMS = [
    "--gamma_d_GB_GluSynapse=101.5",
    "--gamma_p_GB_GluSynapse=199.773931",
    "--a00=1.002", "--a01=2.254638",
    "--a10=1.209857", "--a11=2.396290",
    "--a20=1.127", "--a21=2.500181",
    "--a30=4.214614", "--a31=2.239758",
    "--tau_effca_GB_GluSynapse=278.3177658387",
]


def find_tasks(skip_existing=True):
    tasks, skipped = [], 0
    for pair in sorted(os.listdir(RESULTS_DIR)):
        if not re.fullmatch(r"\d+-\d+", pair):
            continue
        for proto in PROTOCOLS:
            wd = os.path.join(RESULTS_DIR, pair, proto)
            if not os.path.isfile(os.path.join(wd, "prefire_simulation_config.json")):
                skipped += 1
                continue
            out = os.path.join(wd, "simulation_edges_%s_%s.pkl" % (TAG, PARAM_HASH))
            if skip_existing and os.path.isfile(out):
                skipped += 1
                continue
            tasks.append((pair, proto, wd))
    return tasks, skipped


def _worker(task):
    pair, proto, wd = task
    t0 = time.time()
    log_path = os.path.join(wd, "%s_pool_%s_%s.log" % (TAG, proto, pair))
    cmd = ([sys.executable, PAIRRUNNER] + FIT_PARAMS +
           ["--param_hash=%s_%s" % (TAG, PARAM_HASH),
            "--cpre-cpost-cache=%s" % CACHE,
            "--circuit-config=%s" % CIRCUIT])
    try:
        with open(log_path, "w") as logf:
            subprocess.run(cmd, cwd=wd, stdout=logf, stderr=subprocess.STDOUT, check=True)
        ok = os.path.isfile(os.path.join(
            wd, "simulation_edges_%s_%s.pkl" % (TAG, PARAM_HASH)))
        return pair, proto, ok, time.time() - t0, None if ok else "ran but no pkl"
    except subprocess.CalledProcessError as e:
        return pair, proto, False, time.time() - t0, "exit %d" % e.returncode
    except Exception as e:  # noqa: BLE001
        return pair, proto, False, time.time() - t0, str(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=48)
    ap.add_argument("--no-skip-existing", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    tasks, skipped = find_tasks(skip_existing=not args.no_skip_existing)
    print("emodel  : %s" % CIRCUIT)
    print("cache   : %s (exists=%s)" % (CACHE, os.path.exists(CACHE)))
    print("tasks   : %d  (skipped %d)" % (len(tasks), skipped))
    print("workers : %d" % args.workers)
    if args.dry_run:
        for t in tasks[:6]:
            print("  [dry-run] %s / %s" % (t[0], t[1]))
        print("  ... %d total" % len(tasks))
        return
    if not tasks:
        print("Nothing to do.")
        return

    t_start = time.time()
    done = failed = 0
    failures = []
    with multiprocessing.Pool(args.workers, maxtasksperchild=1) as pool:
        for pair, proto, ok, secs, err in pool.imap_unordered(_worker, tasks):
            if ok:
                done += 1
            else:
                failed += 1
                failures.append((pair, proto, err))
            n = done + failed
            rate = n / max(1e-9, time.time() - t_start)
            eta = (len(tasks) - n) / rate if rate else 0
            print("[%d/%d] %s/%s %s (%.0fs)  eta %.1f min" % (
                n, len(tasks), pair, proto,
                "ok" if ok else "FAIL: %s" % err, secs, eta / 60), flush=True)

    print("=" * 70)
    print("Done: %d ok, %d failed in %s" % (
        done, failed, time.strftime("%H:%M:%S", time.gmtime(time.time() - t_start))))
    for pair, proto, err in failures[:30]:
        print("  %s/%s: %s" % (pair, proto, err))
    if len(failures) > 30:
        print("  ... and %d more" % (len(failures) - 30))


if __name__ == "__main__":
    main()
