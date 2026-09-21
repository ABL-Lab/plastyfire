"""
Can apical synapses alone drive LTP?  53 mixed pairs, basal synapses silenced.

Setup
-----
* pairs      : the 53 'mixed' pairs from apical_basal_map.csv (both apical and
               basal synapses present, so silencing basal leaves a real apical
               population -- 142 apical synapses total, median 2 per pair)
* silencing  : PLASTYFIRE_ONLY_LOC=apical, which pushes basal theta_d/theta_p to
               THETA_SILENCE (1e9) so their WATCH never arms and rho_GB' keeps
               them at their initial attractor. NOT theta=-1, which would fire the
               watch immediately and pin dep_GB=pot_GB=1.
* a-params   : all 1.0, per the request for the bare minimum.

CAVEAT ON a=1 -- with every coefficient at 1,
    theta_d = 1*c_pre + 1*c_post = theta_p
so the depression and potentiation thresholds are IDENTICAL. There is no
potentiation-only band: calcium either sits below both (no change) or above both,
where GluSynapse.mod's rho ODE applies `pot_GB*gamma_p*(1-rho) -
dep_GB*(1-pot_GB)*gamma_d*rho`. The (1-pot_GB) factor means potentiation wins
whenever both are active, so a=1 is effectively "potentiation whenever calcium
clears threshold". That is a legitimate question to ask, but it is not a neutral
baseline -- read the result as "can apical calcium clear its own threshold",
not as a balanced LTP/LTD test.

Output: simulation_edges_apicalonly_<hash>.pkl per workdir.
"""

import argparse
import multiprocessing
import os
import subprocess
import sys
import time

import pandas as pd

PLASTYFIRE = "/lustre06/project/6077694/dhuruva/plastyfire"
RESULTS_DIR = os.path.join(
    PLASTYFIRE,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations")
PAIRRUNNER = os.path.join(PLASTYFIRE, "plastyfire/pairrunner_edges_fit.py")
CACHE = os.path.join(PLASTYFIRE, "cpre_cpost_cache/ion_channels_tau278.pkl")
CIRCUIT = os.path.join(PLASTYFIRE,
                       "data/dhuruva_modified_ion_channels_circuit_config.json")
MAP_CSV = os.path.join(PLASTYFIRE, "apical_basal_map.csv")
TAG = "apicalonly"
HASH = "0ce64fa83b85"
PROTOCOLS = ["10Hz_10ms", "10Hz_-10ms"]

# All a-params at 1.0 (the "bare minimum" request). gamma/tau kept at gen-4 values
# so only the threshold coefficients are neutralised.
FIT_PARAMS = [
    "--gamma_d_GB_GluSynapse=101.5",
    "--gamma_p_GB_GluSynapse=199.773931",
    "--a00=1.0", "--a01=1.0",
    "--a10=1.0", "--a11=1.0",
    "--a20=1.0", "--a21=1.0",
    "--a30=1.0", "--a31=1.0",
    "--tau_effca_GB_GluSynapse=278.3177658387",
]


def find_tasks(skip_existing=True):
    df = pd.read_csv(MAP_CSV)
    mixed = df[df.category == "mixed"]
    tasks, skipped = [], 0
    for r in mixed.itertuples():
        pair = "%d-%d" % (r.pregid, r.postgid)
        for proto in PROTOCOLS:
            wd = os.path.join(RESULTS_DIR, pair, proto)
            if not os.path.isfile(os.path.join(wd, "prefire_simulation_config.json")):
                skipped += 1
                continue
            out = os.path.join(wd, "simulation_edges_%s_%s.pkl" % (TAG, HASH))
            if skip_existing and os.path.isfile(out):
                skipped += 1
                continue
            tasks.append((pair, proto, wd, int(r.n_apical), int(r.n_basal)))
    return tasks, skipped


def _worker(task):
    pair, proto, wd, n_ap, n_ba = task
    t0 = time.time()
    log_path = os.path.join(wd, "%s_pool_%s_%s.log" % (TAG, proto, pair))
    env = dict(os.environ)
    env["PLASTYFIRE_ONLY_LOC"] = "apical"      # <- the silencing switch
    cmd = ([sys.executable, PAIRRUNNER] + FIT_PARAMS +
           ["--param_hash=%s_%s" % (TAG, HASH),
            "--cpre-cpost-cache=%s" % CACHE,
            "--circuit-config=%s" % CIRCUIT])
    try:
        with open(log_path, "w") as logf:
            subprocess.run(cmd, cwd=wd, stdout=logf, stderr=subprocess.STDOUT,
                           check=True, env=env)
        ok = os.path.isfile(os.path.join(wd, "simulation_edges_%s_%s.pkl" % (TAG, HASH)))
        return pair, proto, n_ap, n_ba, ok, time.time() - t0, None if ok else "no pkl"
    except subprocess.CalledProcessError as e:
        return pair, proto, n_ap, n_ba, False, time.time() - t0, "exit %d" % e.returncode
    except Exception as e:  # noqa: BLE001
        return pair, proto, n_ap, n_ba, False, time.time() - t0, str(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=48)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-skip-existing", action="store_true")
    args = ap.parse_args()

    tasks, skipped = find_tasks(skip_existing=not args.no_skip_existing)
    print("mixed pairs -> tasks : %d  (skipped %d)" % (len(tasks), skipped))
    print("silencing            : PLASTYFIRE_ONLY_LOC=apical (basal frozen)")
    print("a-params             : all 1.0  => theta_d == theta_p (see docstring)")
    if args.dry_run:
        for t in tasks[:6]:
            print("  [dry-run] %s / %s  (apical=%d basal=%d)" % (t[0], t[1], t[3], t[4]))
        print("  ... %d total" % len(tasks))
        return
    if not tasks:
        print("Nothing to do.")
        return

    t0 = time.time()
    done = failed = 0
    with multiprocessing.Pool(args.workers, maxtasksperchild=1) as pool:
        for pair, proto, n_ap, n_ba, ok, secs, err in pool.imap_unordered(_worker, tasks):
            done, failed = (done + 1, failed) if ok else (done, failed + 1)
            n = done + failed
            print("[%d/%d] %s/%s ap=%d ba=%d %s (%.0fs)" % (
                n, len(tasks), pair, proto, n_ap, n_ba,
                "ok" if ok else "FAIL: %s" % err, secs), flush=True)
    print("=" * 70)
    print("Done: %d ok, %d failed in %s" % (
        done, failed, time.strftime("%H:%M:%S", time.gmtime(time.time() - t0))))


if __name__ == "__main__":
    main()
