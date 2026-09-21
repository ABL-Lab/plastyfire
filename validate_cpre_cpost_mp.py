#!/usr/bin/env python3
"""
Validate c_pre/c_post against the published Ecker 2023 thresholds.

WHY THIS EXISTS. There are three implementations of "compute c_pre/c_post" in this
repo and they do not agree:

  1. plastyfire/simulator.py    c_pre_finder / c_post_finder  — the CANONICAL pair,
     driven by thresholdfinder.py; this lineage produced the published thresholds.
  2. compute_thresholds.py      inline reimplementation; scores median 0.97-0.98 /
     ~89% within 20% against the published circuit (757 synapses).
  3. simulator_edges.py         _find_cpre_cpost — built every cpre_cpost_cache/*.pkl.
     Scores ~0.97-1.01 here; the Jul 6 cache scores 16% on apicals (median 0.49).

--engine selects which one runs, so they can be compared on identical inputs against
the same external reference. That is the only way to tell an implementation bug from
an input difference.

KNOWN DIFFERENCES between (1) and (3), found by reading them side by side:

  a. c_post stimulus. (1) uses h.TStim + tstim.train(offset, duration, amp, freq,
     width) with the stimulus dict from thresholdfinder.py. (3) hardcodes a single
     h.IClamp with dur=3.0 and scrapes amp from the workdir config's amp_start.
     thresholdfinder.py searches amp and width TOGETHER (pulse_width in [1.5, 3, 5],
     spike_threshold_finder returns the amp that fires exactly 1 AP at that width),
     so an amp is only valid with its own width. Worse, the stored single_cells pkl
     holds nspikes=5 / freq=10 — the 5-pulse induction train — so that amp was never
     a single-spike threshold at all.
  b. c_pre synapse setup. (1) and (2) set rho0_GB=1, Use_p=1, Use=Use_p,
     gmax0_AMPA=gmax_p_AMPA. (3) additionally sets Use_GB=1 and gmax_AMPA — neither
     appears in either reference implementation.

Comparison is on theta because the reference stores only theta; c_pre/c_post are then
recovered by inverting the a-param matrix (as DEES_cell_packages/tests/test_CPre.py
does), so a c_pre error and a c_post error are separable rather than hidden inside a
weighted sum. For apicals theta_p is c_pre-dominated (a30=5.24 vs a31=1.78) while
theta_d carries far more c_post weight (a21=2.46 vs a20=1.13).

Usage:
    # canonical engine vs published thresholds — the reference check
    python validate_cpre_cpost_mp.py --engine canonical --limit 5 --workers 5

    # the edges engine that built the caches, same inputs
    python validate_cpre_cpost_mp.py --engine edges --limit 5 --workers 5

    # reproduce a cache (same-code check, not a correctness check)
    python validate_cpre_cpost_mp.py --reference cache --engine edges --workers 60 \
        --cache cpre_cpost_cache/ion_channels_tau278.pkl \
        --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json
"""

import argparse
import json
import os
import pickle
import sys
import time

import h5py
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from precompute_cpre_cpost import (  # noqa: E402
    PARAM_PRESETS, find_workdirs, _NoDaemonPool,
)

EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
NODE_POP = "S1nonbarrel_neurons"
RESULTS_DIR = "/project/ctb-emuller/dhuruva/plastyfire/refitting_results"
APICAL_SECTION_TYPE = 3

# Published Ecker 2023 circuit — read-only, owned by emuller, untouched by this
# project. Its circuit_config.json names the same edges.h5 the theta is read from, so
# synapse properties and reference thresholds are one self-consistent source.
ECKER = "/project/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker"
OLD_EDGES = os.path.join(
    ECKER, "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical", "edges.h5")
OLD_CONFIG = os.path.join(ECKER, "circuit_config.json")
SINGLE_CELLS = os.path.join(
    HERE, "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP",
    "simulations", "single_cells")

# Full-precision a-params from config/hexO1_v7.yaml via test_CPre.py — the values the
# published theta was built with. compute_thresholds.py carries these rounded to 4
# significant figures, a ~0.02% error: irrelevant at 20% tolerance, not at 0.001.
CHINDEMI_A = {"a00": 1.0018099627, "a01": 1.9535568661,
              "a10": 1.1593870631, "a11": 2.4827933785,
              "a20": 1.1267343377, "a21": 2.4559713296,
              "a30": 5.2356637311, "a31": 1.7822214534}

PROTOS = ["10Hz_-10ms", "10Hz_10ms", "10Hz_5ms", "10Hz_-30ms",
          "10Hz_30ms", "10Hz_-50ms", "10Hz_50ms"]


# ── reference maths ─────────────────────────────────────────────────────────────

def theta_from_c(c_pre, c_post, apical, a=CHINDEMI_A):
    if apical:
        return a["a20"] * c_pre + a["a21"] * c_post, a["a30"] * c_pre + a["a31"] * c_post
    return a["a00"] * c_pre + a["a01"] * c_post, a["a10"] * c_pre + a["a11"] * c_post


def _inv_matrices(a=CHINDEMI_A):
    """Inverses of the 2x2 a-param matrices — same as test_CPre.py's inv_params_fn."""
    return (np.linalg.inv([[a["a00"], a["a01"]], [a["a10"], a["a11"]]]),
            np.linalg.inv([[a["a20"], a["a21"]], [a["a30"], a["a31"]]]))


def c_from_theta(theta_d, theta_p, apical, inv):
    """Invert (theta_d, theta_p) -> (c_pre, c_post) to separate the two terms."""
    m = inv[1] if apical else inv[0]
    return m[0, 0] * theta_d + m[0, 1] * theta_p, m[1, 0] * theta_d + m[1, 1] * theta_p


def read_theta(edges_h5, syn_ids):
    """Batch-read theta_d/theta_p/section_type. One sorted read for every synapse."""
    ids = np.unique(np.asarray(syn_ids, dtype=np.int64))
    p = f"edges/{EDGE_POP}/0"
    with h5py.File(edges_h5, "r") as f:
        td, tp = f[f"{p}/theta_d"][ids], f[f"{p}/theta_p"][ids]
        st = f[f"{p}/afferent_section_type"][ids]
    return {int(s): (float(a), float(b), int(c))
            for s, a, b, c in zip(ids, td, tp, st)}


def load_stimulus(post_gid):
    """The stimulus dict c_post_finder expects, from the single-cell pkl.

    Carries width/amp as the matched pair thresholdfinder.py searched for, instead of
    pairing a scraped amp with a hardcoded 3.0 ms duration.
    """
    pkl = os.path.join(SINGLE_CELLS, f"{post_gid}.pkl")
    if not os.path.exists(pkl):
        return None
    with open(pkl, "rb") as f:
        data = pickle.load(f)
    d = data[list(data.keys())[0]]
    return {k: d[k] for k in ("nspikes", "freq", "width", "offset", "amp")}


# ── engines ─────────────────────────────────────────────────────────────────────

def _worker(args):
    """Run one pair through the selected engine. Top-level so it stays picklable."""
    workdir, pre_gid, post_gid, fit_params, circuit_config, fixhp, engine = args
    try:
        if engine == "edges":
            import plastyfire.simulator_edges as sim_mod
            res = sim_mod.compute_cpre_cpost_for_workdir(
                workdir, fit_params=fit_params, node_pop=NODE_POP, edge_pop=EDGE_POP,
                fixhp=fixhp, circuit_config=circuit_config)
            if res is None:
                return pre_gid, post_gid, None, None, "returned None"
            n_sp = res.get("n_post_spikes")
            # A c_post sim in which the cell never fired measures nothing: c_post
            # falls back to the resting effcai_GB floor (~3e-5 vs ~0.07). Surface it
            # as an explicit failure instead of letting it pass as a small number.
            if n_sp == 0:
                return (pre_gid, post_gid, res["c_pre"], res["c_post"],
                        f"POST CELL DID NOT SPIKE at amp={res.get('pulse_amp')} nA "
                        f"— c_post invalid")
            return pre_gid, post_gid, res["c_pre"], res["c_post"], None

        # canonical: plastyfire/simulator.py's own finders, driven exactly as
        # thresholdfinder.py drives them.
        from plastyfire.simulator import c_pre_finder, c_post_finder
        stimulus = load_stimulus(post_gid)
        if stimulus is None:
            return pre_gid, post_gid, None, None, f"no single_cells pkl for {post_gid}"
        # c_pre_finder takes a SIMULATION config (bluecellulab needs its "run"
        # section), not a circuit config. Point the workdir's sim config at the
        # requested circuit instead of passing the circuit config itself.
        sim_config = os.path.join(workdir, "prefire_simulation_config.json")
        if circuit_config:
            with open(sim_config) as f:
                cfg = json.load(f)
            cfg["network"] = os.path.abspath(circuit_config)
            ns = os.path.join(workdir, "node_sets.json")
            if os.path.exists(ns):
                cfg["node_sets_file"] = os.path.abspath(ns)
            cfg.pop("reports", None)
            import tempfile
            tmp = tempfile.NamedTemporaryFile(
                mode="w", suffix=".json", delete=False,
                dir=os.environ.get("SLURM_TMPDIR", tempfile.gettempdir()))
            json.dump(cfg, tmp)
            tmp.close()
            sim_config = tmp.name
        c_pre = c_pre_finder(sim_config, fit_params, None, pre_gid, post_gid,
                             node_pop=NODE_POP, edge_pop=EDGE_POP, fixhp=fixhp)
        c_post = c_post_finder(sim_config, fit_params, None, pre_gid, post_gid,
                               stimulus, node_pop=NODE_POP, edge_pop=EDGE_POP,
                               fixhp=fixhp)
        return pre_gid, post_gid, c_pre, c_post, None
    except Exception as e:  # noqa: BLE001 — one bad pair must not kill the sweep
        return pre_gid, post_gid, None, None, f"{type(e).__name__}: {e}"


# ── reporting ───────────────────────────────────────────────────────────────────

def _stats(g, col):
    r = g[col].replace([np.inf, -np.inf], np.nan).dropna()
    if r.empty:
        return None
    dev = (r - 1).abs()
    return len(r), r.median(), dev.quantile(0.95), dev.max()


def _table(df, pairs, title, cols):
    print(f"\n{title}")
    print(f"  {'loc':>7s} {'n':>5s} | " +
          " | ".join(f"{a + ' med':>12s} {'p95dev':>7s} {'maxdev':>7s}" for a, _ in cols))
    for loc, g in df.groupby("loc"):
        parts = []
        for _, c in cols:
            s = _stats(g, c)
            parts.append("(none)".rjust(28) if s is None
                         else f"{s[1]:>12.4f} {s[2]:>7.3f} {s[3]:>7.3f}")
        print(f"  {loc:>7s} {len(g):>5d} | " + " | ".join(parts))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", choices=["canonical", "edges"], default="canonical",
                    help="canonical: simulator.py's c_pre_finder/c_post_finder (the "
                         "lineage behind the published theta). edges: "
                         "simulator_edges._find_cpre_cpost (built the caches).")
    ap.add_argument("--reference", choices=["published", "cache"], default="published",
                    help="published: compare theta against the Ecker circuit. "
                         "cache: reproduce a cpre/cpost pkl instead.")
    ap.add_argument("--cache", help="cache pkl (required for --reference cache)")
    ap.add_argument("--edges-h5", default=OLD_EDGES)
    ap.add_argument("--circuit-config", default=None,
                    help="default: the published Ecker circuit_config.json")
    ap.add_argument("--params", default="chindemi_aparams", choices=list(PARAM_PRESETS))
    ap.add_argument("--tau-only", action="store_true", default=True,
                    help="apply only tau_effca, as compute_thresholds.py does (default)")
    ap.add_argument("--full-params", dest="tau_only", action="store_false",
                    help="apply the whole preset (gamma_d/gamma_p too)")
    ap.add_argument("--results-dir", default=RESULTS_DIR)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-fixhp", action="store_true")
    ap.add_argument("--tol", type=float, default=0.20)
    ap.add_argument("--csv", default=None)
    args = ap.parse_args()

    if args.reference == "cache" and not args.cache:
        ap.error("--reference cache requires --cache")

    cache = None
    if args.cache:
        with open(args.cache, "rb") as f:
            cache = pickle.load(f)

    if args.circuit_config is None and args.reference == "published":
        args.circuit_config = OLD_CONFIG

    fit_params = dict(PARAM_PRESETS[args.params])
    if args.tau_only:
        fit_params = {"tau_effca_GB_GluSynapse": fit_params["tau_effca_GB_GluSynapse"]}
    fixhp = not args.no_fixhp

    pairs = {}
    for wd in find_workdirs(args.results_dir):
        parts = os.path.basename(os.path.dirname(wd)).split("-")
        if len(parts) != 2:
            continue
        try:
            key = (int(parts[0]), int(parts[1]))
        except ValueError:
            continue
        if cache is not None and key not in cache:
            continue
        pairs.setdefault(key, wd)

    todo = sorted(pairs.items())[:args.limit] if args.limit else sorted(pairs.items())
    if not todo:
        print(f"No pairs found under {args.results_dir}")
        return 1

    print(f"engine    : {args.engine}")
    print(f"reference : {args.reference}"
          + (f" ({os.path.basename(args.edges_h5)})" if args.reference == "published"
             else f" ({args.cache})"))
    print(f"config    : {args.circuit_config or '(workdir prefire config)'}")
    print(f"params    : {args.params}" + (" -> tau_effca only" if args.tau_only else ""))
    print(f"pairs     : {len(todo)}   workers: {args.workers}   fixhp: {fixhp}")
    if args.reference == "published":
        print(f"a-params  : full-precision hexO1_v7 (a30={CHINDEMI_A['a30']})")
    print("=" * 78)

    work = [(wd, k[0], k[1], fit_params, args.circuit_config, fixhp, args.engine)
            for k, wd in todo]

    t0 = time.time()
    got, failures = {}, []
    with _NoDaemonPool(processes=min(args.workers, len(work))) as pool:
        for i, (pre, post, c_pre, c_post, err) in enumerate(
                pool.imap_unordered(_worker, work), 1):
            if err:
                failures.append(((pre, post), err))
                print(f"  [{i}/{len(work)}] ({pre}, {post})  FAILED: {err}", flush=True)
            else:
                got[(pre, post)] = {"c_pre": c_pre, "c_post": c_post}
                print(f"  [{i}/{len(work)}] ({pre}, {post})  {len(c_pre)} synapses",
                      flush=True)

    if not got:
        print("\nEvery pair failed.")
        for k, e in failures[:10]:
            print(f"  {k}: {e}")
        return 1

    all_syn = [int(s) for v in got.values() for s in v["c_pre"]]
    rows = []

    if args.reference == "published":
        ref = read_theta(args.edges_h5, all_syn)
        inv = _inv_matrices()
        for key, new in got.items():
            for syn in new["c_pre"]:
                syn = int(syn)
                td_ref, tp_ref, st = ref[syn]
                apical = st == APICAL_SECTION_TYPE
                cp, cq = float(new["c_pre"][syn]), float(new["c_post"][syn])
                td_new, tp_new = theta_from_c(cp, cq, apical)
                cp_ref, cq_ref = c_from_theta(td_ref, tp_ref, apical, inv)
                rows.append({
                    "pair": f"{key[0]}-{key[1]}", "syn_id": syn,
                    "loc": "apical" if apical else "basal",
                    "c_pre_new": cp, "c_pre_ref": cp_ref,
                    "c_post_new": cq, "c_post_ref": cq_ref,
                    "cpre_ratio": cp / cp_ref if cp_ref else np.nan,
                    "cpost_ratio": cq / cq_ref if cq_ref else np.nan,
                    "theta_d_new": td_new, "theta_d_ref": td_ref,
                    "theta_p_new": tp_new, "theta_p_ref": tp_ref,
                    "td_ratio": td_new / td_ref if td_ref > 0 else np.nan,
                    "tp_ratio": tp_new / tp_ref if tp_ref > 0 else np.nan,
                })
        pri = ("td_ratio", "tp_ratio")
    else:
        with h5py.File(args.edges_h5, "r") as f:
            ids = np.unique(np.asarray(all_syn, dtype=np.int64))
            st = f[f"edges/{EDGE_POP}/0/afferent_section_type"][ids]
        sec = dict(zip(ids.tolist(), st.tolist()))
        for key, new in got.items():
            old = cache.get(key)
            if old is None:
                continue
            for syn in new["c_pre"]:
                syn = int(syn)
                if syn not in old["c_pre"]:
                    continue
                cp, cq = float(new["c_pre"][syn]), float(new["c_post"][syn])
                cp_o, cq_o = float(old["c_pre"][syn]), float(old["c_post"][syn])
                rows.append({
                    "pair": f"{key[0]}-{key[1]}", "syn_id": syn,
                    "loc": "apical" if sec[syn] == APICAL_SECTION_TYPE else "basal",
                    "c_pre_new": cp, "c_pre_ref": cp_o,
                    "c_post_new": cq, "c_post_ref": cq_o,
                    "cpre_ratio": cp / cp_o if cp_o else np.nan,
                    "cpost_ratio": cq / cq_o if cq_o else np.nan,
                })
        pri = ("cpre_ratio", "cpost_ratio")

    df = pd.DataFrame(rows)
    if df.empty:
        print("\nNothing to compare.")
        return 1
    df["ok"] = ((df[pri[0]] - 1).abs() <= args.tol) & \
               ((df[pri[1]] - 1).abs() <= args.tol)

    print("\n" + "=" * 78)
    print(f"{len(df)} synapses / {len(got)} pairs in {time.time() - t0:.0f}s "
          f"({len(failures)} failed)")

    if args.reference == "published":
        _table(df, got, "theta, recomputed / published:",
               [("theta_d", "td_ratio"), ("theta_p", "tp_ratio")])
    _table(df, got, "c_pre / c_post, recomputed / reference:",
           [("c_pre", "cpre_ratio"), ("c_post", "cpost_ratio")])

    print("\nabsolute medians:")
    for loc, g in df.groupby("loc"):
        print(f"  {loc:>7s}  c_pre  new={g['c_pre_new'].median():.6f} "
              f"ref={g['c_pre_ref'].median():.6f}  |  "
              f"c_post new={g['c_post_new'].median():.8f} "
              f"ref={g['c_post_ref'].median():.8f}")

    n_pass = int(df["ok"].sum())
    print(f"\nWithin +/-{args.tol:.0%}: {n_pass}/{len(df)} ({n_pass / len(df):.1%})")

    worst = df.reindex((df[pri[0]] - 1).abs().sort_values(ascending=False).index)
    print(f"\n10 worst {pri[0]}:")
    print(worst[["pair", "syn_id", "loc", "c_pre_new", "c_pre_ref",
                 "c_post_new", "c_post_ref", pri[0]]].head(10).to_string(index=False))

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"\nWrote {args.csv}")
    if failures:
        print(f"\n{len(failures)} failed:")
        for k, e in failures[:10]:
            print(f"  {k}: {e}")

    print(f"\n{'MATCH' if n_pass == len(df) else 'MISMATCH'}")
    return 0 if n_pass == len(df) else 1


if __name__ == "__main__":
    sys.exit(main())
