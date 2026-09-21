"""
Extracellular calcium sweep of the STDP induction sims, DE fit #2 params.

Runs the same 700 pair x protocol bluecellulab simulations as run_de_fit2_pool.py
at several extracellular Ca levels (default 3.0 / 1.8 / 1.3 mM vs the 2.0 mM the
model was fitted at), recording rho and nothing else.

Thresholds are rebuilt at every calcium level
---------------------------------------------
theta_d = a00*c_pre + a01*c_post (basal; a20/a21 apical) and theta_p likewise,
and c_pre/c_post are themselves calcium-dependent: cao_CR enters GluSynapse.mod
as Pf_NMDA = (4*cao_CR)/(4*cao_CR + (1/1.38)*120)*0.6 (line 323) and through
nernst(cai_CR, cao_CR, 2) in Eca_syn (line 328). Pf_NMDA alone moves -33% at
1.3 mM and +44% at 3.0 mM relative to 2.0 mM. Reusing the 2.0 mM cache would
therefore hold the thresholds fixed while the calcium transients moved under
them, making every synapse spuriously easy to potentiate at high Ca and hard at
low Ca -- measuring the threshold error, not the calcium biology.

So each Ca level gets its OWN cpre/cpost cache, built by precompute_cpre_cpost.py
at that Ca, and the a-params (which are Ca-independent fitted constants) turn it
into that level's thresholds at runtime. Only tau_effca_GB affects effcai_GB
(mod line 360), and the `chindemi` preset already carries fit #2's exact
tau=278.3177658387, so it is the correct preset for the cache mini-sims.

edges.h5 is never touched: FIT_ARGS carries a-params, so the BCL path computes
theta at runtime from the cache (_apply_theta_from_a_params) rather than reading
edges.h5. The concurrently running neurodamus jobs are therefore unaffected.

The EPSP basis IS reused across levels: the reported quantity is a ratio
(after / before) at one Ca, so the Ca-dependent release scaling is a common
factor that divides out.

Output (one small pkl per pair x protocol per Ca, nothing written into the
workdirs except the usual edges_prefire.log):

    ca_scan/ca_1.30mM/183120-189281__10Hz_-10ms.pkl
        {global_ids, initial_rho, final_rho, n_post_spikes, n_post_spikes_exp}

~1 kB each instead of 520 kB, and --lean also suppresses the 275 MB
simulation_traces.pkl and the per-workdir rho.h5 that the 2.0 mM run left behind.

Usage (inside an salloc/sbatch allocation) -- caches first, then sims:
    python run_ca_scan_pool.py --workers 60 --build-caches
    python run_ca_scan_pool.py --workers 60
    python run_ca_scan_pool.py --workers 60 --dry-run
"""

import argparse
import multiprocessing
import os
import pickle
import shutil
import subprocess
import sys
import tempfile
import time

# Single source of truth for the fitted params: DE fit #2.
from run_de_fit2_pool import FIT_ARGS, PAIRRUNNER, PROTOCOLS, RESULTS_DIR

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT = os.path.join(PLASTYFIRE_ROOT, "ca_scan")
CA_LEVELS = [3.0, 1.8, 1.3]

# --param_hash only names the default output file, which --out replaces.
# --cpre-cpost-cache is dropped here and re-added per Ca level by _worker: the
# 2.0 mM cache in FIT_ARGS is invalid at any other calcium.
BASE_ARGS = [a for a in FIT_ARGS
             if not a.startswith(("--param_hash", "--cpre-cpost-cache"))]

# tau_effca is the only HOC global that affects effcai_GB (GluSynapse.mod:360);
# gamma_d/gamma_p act on rho, which never moves in the mini-sims (theta=-1).
# `chindemi` carries tau=278.3177658387 == fit #2's tau.
# Peak RSS per bluecellulab induction process. run_de_fit2_pool's own note says
# "keep workers * 4 GB under the job's --mem"; 4.5 adds the headroom that note
# lacked -- at 60 workers x 4 GB = 240 GB against a 250 G allocation the node hit
# FreeMem=583 MB and CPULoad=185 (thrashing, not progress).
MEM_PER_PROC_GB = 4.5
# Leave 20% for page cache and the OS. Filling the cgroup to the brim is what
# produced FreeMem=583 MB: the kernel then evicts the NEURON .so pages it is
# actively executing and re-reads them from lustre, so the node crawls.
MEM_HEADROOM = 0.8

# One thread per process. Each worker is already a whole core; letting numpy /
# OpenBLAS also fan out multiplies runnable threads by the core count and turns
# oversubscription into thrashing.
THREAD_ENV = {v: "1" for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                               "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                               "VECLIB_MAXIMUM_THREADS")}

CACHE_PRESET = "chindemi"
PRECOMPUTE = os.path.join(PLASTYFIRE_ROOT, "precompute_cpre_cpost.py")
CIRCUIT_CONFIG = os.path.join(
    PLASTYFIRE_ROOT, "data/dhuruva_modified_ion_channels_circuit_config.json")


def ca_tag(ca):
    return f"ca_{ca:.2f}mM"


def cache_path(ca):
    """Per-Ca cpre/cpost cache. Separate files: precompute_cpre_cpost.py resumes
    from whatever the output file already holds, so a shared path would silently
    keep the previous level's entries."""
    return os.path.join(PLASTYFIRE_ROOT, "cpre_cpost_cache",
                        f"fit2_{ca_tag(ca)}.pkl")


def build_caches(ca_levels, workers):
    """Run the c_pre/c_post mini-sims once per Ca level (98 usable pairs each)."""
    for ca in ca_levels:
        out = cache_path(ca)
        cmd = [sys.executable, PRECOMPUTE,
               "--params", CACHE_PRESET,
               # find_workdirs globs <root>/fitting/*/seed*/*_STDP/simulations/...,
               # so it needs the root ABOVE fitting/, not the simulations dir.
               "--results-dir", RESULTS_DIR.rstrip("/").split("/fitting/")[0],
               "--output", out,
               "--circuit-config", CIRCUIT_CONFIG,
               f"--extracellular-calcium={ca}",
               "--workers", str(workers)]
        print(f"\n=== cpre/cpost cache @ {ca} mM -> {out} ===", flush=True)
        print("  " + " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)


def find_tasks(ca_levels, skip_existing):
    tasks, skipped = [], 0
    for ca in ca_levels:
        for pair in sorted(os.listdir(RESULTS_DIR)):
            pair_dir = os.path.join(RESULTS_DIR, pair)
            if not os.path.isdir(pair_dir):
                continue
            for proto in PROTOCOLS:
                wd = os.path.join(pair_dir, proto)
                if not os.path.isfile(os.path.join(wd, "prefire_simulation_config.json")):
                    skipped += 1
                    continue
                out = os.path.join(OUT_ROOT, ca_tag(ca), f"{pair}__{proto}.pkl")
                if skip_existing and os.path.isfile(out):
                    skipped += 1
                    continue
                tasks.append((ca, pair, proto, wd, out))
    return tasks, skipped


def _worker(task):
    """One simulation, own process, cwd=workdir. Logs kept only on failure."""
    ca, pair, proto, wd, out = task
    t0 = time.time()
    cmd = [sys.executable, PAIRRUNNER] + BASE_ARGS + [
        f"--cpre-cpost-cache={cache_path(ca)}",
        f"--extracellular-calcium={ca}", "--lean", f"--out={out}"]
    fd, tmp_log = tempfile.mkstemp(suffix=".log")
    os.close(fd)
    try:
        env = dict(os.environ, **THREAD_ENV)
        with open(tmp_log, "w") as logf:
            subprocess.run(cmd, cwd=wd, stdout=logf, stderr=subprocess.STDOUT,
                           check=True, env=env)
        err = None if os.path.isfile(out) else "ran but no pkl written"
    except subprocess.CalledProcessError as e:
        err = f"exit {e.returncode}"
    except Exception as e:                                    # pragma: no cover
        err = str(e)
    if err is None:
        os.unlink(tmp_log)
    else:
        keep = os.path.join(OUT_ROOT, ca_tag(ca), "_failed", f"{pair}__{proto}.log")
        os.makedirs(os.path.dirname(keep), exist_ok=True)
        shutil.move(tmp_log, keep)
        err = f"{err} (see {keep})"
    return ca, pair, proto, err is None, time.time() - t0, err


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=40,
                    help="Parallel simulations (default: 40). ~4.5 GB peak each — "
                         "memory, not core count, is the limit.")
    ap.add_argument("--ca", type=float, nargs="+", default=CA_LEVELS,
                    help=f"Extracellular Ca levels in mM (default: {CA_LEVELS})")
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force-workers", action="store_true",
                    help="Skip the memory guard (expect thrashing).")
    ap.add_argument("--build-caches", action="store_true",
                    help="Build the per-Ca c_pre/c_post caches (mini-sims) and exit. "
                         "Must be run before the induction sims.")
    args = ap.parse_args()

    if args.build_caches:
        build_caches(args.ca, args.workers)
        print("\nCaches built:")
        for ca in args.ca:
            pth = cache_path(ca)
            n = len(pickle.load(open(pth, "rb"))) if os.path.isfile(pth) else 0
            print(f"  {ca_tag(ca)}: {n} pairs  {pth}")
        return

    # Refuse to silently fall back to mini-sims (or to the wrong Ca) mid-sweep.
    missing = [ca for ca in args.ca if not os.path.isfile(cache_path(ca))]
    if missing:
        sys.exit(f"No c_pre/c_post cache for {missing} mM.\n"
                 f"Thresholds are calcium-dependent, so each level needs its own.\n"
                 f"Run first:  python {os.path.basename(__file__)} "
                 f"--workers {args.workers} --build-caches")

    # Memory, not CPU count, is the binding constraint: 63 cores but only enough
    # RAM for ~40 concurrent bluecellulab processes.
    mem_mb = int(os.environ.get("SLURM_MEM_PER_NODE", 0))
    if mem_mb:
        safe = max(1, int(mem_mb / 1024 * MEM_HEADROOM / MEM_PER_PROC_GB))
        if args.workers > safe:
            print(f"!! --workers {args.workers} needs "
                  f"~{args.workers * MEM_PER_PROC_GB:.0f} GB but the job has "
                  f"{mem_mb/1024:.0f} GB.")
            if args.force_workers:
                print("   --force-workers given; continuing (expect thrashing).")
            else:
                sys.exit(f"   Capping is safer: rerun with --workers {safe}.\n"
                         f"   (override with --force-workers)")
        else:
            print(f"memory   : {mem_mb/1024:.0f} GB, "
                  f"{args.workers} x {MEM_PER_PROC_GB} GB = "
                  f"{args.workers * MEM_PER_PROC_GB:.0f} GB — ok")

    tasks, skipped = find_tasks(args.ca, args.skip_existing)
    print(f"ca levels: {args.ca} mM   (thresholds rebuilt per level)")
    for ca in args.ca:
        print(f"  {ca_tag(ca)} cache: {len(pickle.load(open(cache_path(ca),'rb')))} pairs")
    print(f"tasks    : {len(tasks)}  (skipped {skipped})")
    print(f"workers  : {args.workers}")
    print(f"out      : {OUT_ROOT}/")
    if args.dry_run:
        for t in tasks[:5]:
            print(f"  [dry-run] {t[0]} mM  {t[1]}/{t[2]} -> {t[4]}")
        print(f"  … {len(tasks)} total")
        return
    if not tasks:
        print("Nothing to do.")
        return

    for ca in args.ca:
        os.makedirs(os.path.join(OUT_ROOT, ca_tag(ca)), exist_ok=True)

    t_start = time.time()
    done = failed = 0
    failures = []
    with multiprocessing.Pool(args.workers, maxtasksperchild=1) as pool:
        for ca, pair, proto, ok, secs, err in pool.imap_unordered(_worker, tasks):
            done, failed = done + ok, failed + (not ok)
            if not ok:
                failures.append((ca, pair, proto, err))
            n = done + failed
            rate = n / max(1e-9, time.time() - t_start)
            print(f"[{n}/{len(tasks)}] {ca}mM {pair}/{proto} "
                  f"{'ok' if ok else 'FAIL: ' + str(err)} ({secs:.0f}s)  "
                  f"eta {(len(tasks) - n) / rate / 60:.1f} min", flush=True)

    print("=" * 70)
    print(f"Done: {done} ok, {failed} failed in "
          f"{time.strftime('%H:%M:%S', time.gmtime(time.time() - t_start))}")
    for ca in args.ca:
        d = os.path.join(OUT_ROOT, ca_tag(ca))
        n = len([f for f in os.listdir(d) if f.endswith(".pkl")]) if os.path.isdir(d) else 0
        print(f"  {ca_tag(ca)}: {n} pkls")
    for ca, pair, proto, err in failures[:40]:
        print(f"  FAIL {ca}mM {pair}/{proto}: {err}")
    if len(failures) > 40:
        print(f"  … and {len(failures) - 40} more")


if __name__ == "__main__":
    main()
