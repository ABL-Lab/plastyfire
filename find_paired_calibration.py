#!/usr/bin/env python3
"""
Search for a postsynaptic pulse amplitude, WITH the presynaptic synapse
attached, that reliably reproduces the calibrated induction spike count
(nspikes * nreps) for pairs where simwriter.py's isolated-cell calibration
(check_electrical_constraint) undershoots.

Background: simwriter.py calibrates the induction pulse amplitude per POST
cell alone (current clamp, no synapse), once per gid, deliberately -- see the
docstring on check_electrical_constraint. The real pairing simulation adds
synaptic drive from the specific presynaptic partner on top of that pulse,
which can push some pairs a few spikes over the calibrated target. This is
what simulator_edges.py's induction guardrail catches.

Checked against the 4 param sets already run through run_de_fit2_pool.py
(bdbf06f915d0, 21f3bec952eb, 283d79c04294, b8c7ff3ecf0a): 26 of the ~26-28
guardrail failures are IDENTICAL across all four, including DE fit #2, which
has a completely different theta formula and gamma dynamics. So this is a
calibration property of ~19 specific pairs, not an artifact of the
plasticity parameters being tested -- meaning this search needs to run once,
not once per parameter sweep.

Two things ruled out empirically before writing this (see conversation):
  - Checking only one dt (e.g. +5ms) misses ~53% of affected pairs: 9 of 17
    fail ONLY at -10ms or +10ms, never at +5ms. Sign/size of dt controls when
    the presynaptic EPSP lands relative to the post pulse/AHP, and that's
    pair-specific. -> search every dt where THAT pair actually failed.
  - Checking only repetition 1 misses ~46% of the (pair,protocol) failures:
    large excesses (59-60 spikes) show an extra spike every rep from the
    start, but marginal ones (51-53 spikes) are a single stochastic extra
    spike that can land anywhere across the full 10-rep train.
    -> must run the FULL induction, not a truncated proxy.

Given both, this reuses the real, already-validated induction path
(pairrunner_edges_fit.py --force) rather than new simulation code: for each
(pair, protocol) it copies the real workdir's config into a scratch dir,
scales every pulse's amp_start by a candidate fraction, runs the full
induction, and reads n_post_spikes back from the resulting pkl. Binary
search on the fraction, since firing rate is expected to be monotonic in
pulse amplitude (the same assumption spike_threshold_finder itself makes).

Usage:
    python find_paired_calibration.py --dry-run
    python find_paired_calibration.py --pilot 3 --workers 3      # sanity check
    python find_paired_calibration.py --workers 8                # every target
"""
import argparse
import json
import os
import pickle
import shutil
import subprocess
import sys

PLASTYFIRE_ROOT = "/lustre06/project/6077694/dhuruva/plastyfire"
DEFAULT_STDP_NAME = "L5TTPC_L5TTPC_STDP"
# Reference hash to mine for the known-bad (pair, protocol) target list.
# a-params are irrelevant to spike count (proven param-independent -- see
# above), so any completed hash works; must have a completed run to mine from
# for whichever --stdp-name is in play (bdbf06f915d0 only exists for
# L5TTPC_L5TTPC_STDP; use --reference-hash b8c7ff3ecf0a for L23PC_L5TTPC_STDP
# once that run completes).
DEFAULT_REFERENCE_HASH = "bdbf06f915d0"


def sim_root_for(stdp_name):
    return os.path.join(
        PLASTYFIRE_ROOT,
        f"refitting_results/fitting/n100/seed19091997/{stdp_name}/simulations",
    )


SIM_ROOT = sim_root_for(DEFAULT_STDP_NAME)      # overridden in main() via --stdp-name
REFERENCE_HASH = DEFAULT_REFERENCE_HASH          # overridden in main() via --reference-hash
SCRATCH_ROOT = os.path.join(PLASTYFIRE_ROOT, "paired_calibration")
PAIRRUNNER = os.path.join(PLASTYFIRE_ROOT, "plastyfire/pairrunner_edges_fit.py")
CACHE = os.path.join(PLASTYFIRE_ROOT, "cpre_cpost_cache/ion_channels_tau278.pkl")
CIRCUIT_CONFIG = os.path.join(
    PLASTYFIRE_ROOT, "data/dhuruva_modified_ion_channels_circuit_config.json"
)

# Spike count is a-param-independent FAR from the firing threshold (proven
# across 4 hashes at the original, too-strong amplitude -- see module
# docstring). It is NOT independent right AT a precisely-tuned threshold: of
# 27 amplitudes found under the "cpreonly" preset below, 8 still missed the
# target under DE fit #2's very different gammas (77.76/299.91) and nonzero
# c_post term once patched into the real workdirs and run for real -- a 1-2
# spike difference invisible while the amplitude was still far too high.
# --fit-preset picks which set drives the search; use whichever hash you are
# actually about to run, not just the original validation preset.
FIT_PRESETS = {
    "cpreonly": [   # bdbf06f915d0 -- the original validation preset (default)
        "--gamma_d_GB_GluSynapse=84.3845", "--gamma_p_GB_GluSynapse=186.5408",
        "--a00=1.004709", "--a01=0.0", "--a10=1.977768", "--a11=0.0",
        "--a20=1.004709", "--a21=0.0", "--a30=1.977768", "--a31=0.0",
    ],
    "defit2": [     # b8c7ff3ecf0a -- DE fit #2's actual parameters
        "--gamma_d_GB_GluSynapse=77.7558", "--gamma_p_GB_GluSynapse=299.9121",
        "--a00=1.003498", "--a01=2.902478", "--a10=1.644558", "--a11=2.764812",
        "--a20=1.003498", "--a21=2.902478", "--a30=1.644558", "--a31=2.764812",
    ],
}

FIT_ARGS = FIT_PRESETS["cpreonly"] + [
    "--tau_effca_GB_GluSynapse=278.3177658387",
    "--param_hash=calib_scan",
    f"--cpre-cpost-cache={CACHE}",
    f"--circuit-config={CIRCUIT_CONFIG}",
    "--force", "--lean",
]

FRAC_LO, FRAC_HI = 0.5, 1.0   # search only downward: every known failure is an EXCESS
MAX_ITERS = 7                  # 0.5/2^7 ~ 0.004 resolution in fraction space


def find_targets():
    """(pair, protocol) tuples with guardrail_forced=True under REFERENCE_HASH."""
    targets = []
    for pair in sorted(os.listdir(SIM_ROOT)):
        pair_dir = os.path.join(SIM_ROOT, pair)
        if not os.path.isdir(pair_dir):
            continue
        for proto in sorted(os.listdir(pair_dir)):
            pkl = os.path.join(pair_dir, proto, f"simulation_edges_{REFERENCE_HASH}.pkl")
            if not os.path.isfile(pkl):
                continue
            with open(pkl, "rb") as f:
                d = pickle.load(f)
            if d.get("guardrail_forced"):
                targets.append((pair, proto, d["n_post_spikes"], d["n_post_spikes_exp"]))
    return targets


def _run_one_amplitude(pair, proto, frac, trial_dir):
    """Copy the real workdir into trial_dir with pulses scaled by frac, run the
    full induction, return n_post_spikes (or None on failure)."""
    src = os.path.join(SIM_ROOT, pair, proto)
    os.makedirs(trial_dir, exist_ok=True)

    with open(os.path.join(src, "prefire_simulation_config.json")) as f:
        cfg = json.load(f)
    pulse_keys = sorted(k for k in cfg["inputs"] if k.startswith("pulse"))
    orig_amps = {k: cfg["inputs"][k]["amp_start"] for k in pulse_keys}
    assert len(set(orig_amps.values())) == 1, \
        f"{pair}/{proto}: pulses do not share one amplitude, unexpected: {orig_amps}"
    for k in pulse_keys:
        cfg["inputs"][k]["amp_start"] = orig_amps[k] * frac
    with open(os.path.join(trial_dir, "prefire_simulation_config.json"), "w") as f:
        json.dump(cfg, f)

    # prefire_prespikes.h5 / node_sets.json are read-only and dt/pair-specific
    # but amplitude-independent -- copy once, never modified.
    shutil.copy(os.path.join(src, "prefire_prespikes.h5"), trial_dir)
    ns = os.path.join(src, "node_sets.json")
    if os.path.isfile(ns):
        shutil.copy(ns, trial_dir)

    out_pkl = os.path.join(trial_dir, "result.pkl")
    log_path = os.path.join(trial_dir, "run.log")
    cmd = [sys.executable, PAIRRUNNER] + FIT_ARGS + ["--out", out_pkl]
    with open(log_path, "w") as logf:
        subprocess.run(cmd, cwd=trial_dir, stdout=logf, stderr=subprocess.STDOUT)

    if not os.path.isfile(out_pkl):
        return None
    with open(out_pkl, "rb") as f:
        return pickle.load(f)["n_post_spikes"]


def apply_patches(results):
    """Write each converged found_frac into the REAL workdir's
    prefire_simulation_config.json, so every future fit (any new hash) gets
    the corrected amplitude automatically instead of rediscovering these same
    failures. The untouched original is preserved alongside it so this is
    reversible; re-running --apply is idempotent (always scales from that
    preserved original, never from an already-patched file, so repeat runs
    can't compound the scaling)."""
    applied = []
    for r in results:
        if not r["converged"]:
            print(f"SKIP (not converged): {r['pair']}/{r['proto']}")
            continue
        wd = os.path.join(SIM_ROOT, r["pair"], r["proto"])
        cfg_path = os.path.join(wd, "prefire_simulation_config.json")
        backup_path = cfg_path + ".orig_amp_backup"

        if os.path.isfile(backup_path):
            with open(backup_path) as f:
                cfg = json.load(f)
            print(f"  {r['pair']}/{r['proto']}: already patched once, "
                  f"re-deriving from preserved original")
        else:
            with open(cfg_path) as f:
                cfg = json.load(f)
            shutil.copy(cfg_path, backup_path)

        pulse_keys = sorted(k for k in cfg["inputs"] if k.startswith("pulse"))
        orig_amps = {k: cfg["inputs"][k]["amp_start"] for k in pulse_keys}
        assert len(set(orig_amps.values())) == 1, \
            f"{r['pair']}/{r['proto']}: pulses do not share one amplitude: {orig_amps}"
        orig_amp = next(iter(orig_amps.values()))
        new_amp = orig_amp * r["found_frac"]
        for k in pulse_keys:
            cfg["inputs"][k]["amp_start"] = new_amp
        with open(cfg_path, "w") as f:
            json.dump(cfg, f)

        applied.append(dict(pair=r["pair"], proto=r["proto"], orig_amp=orig_amp,
                            new_amp=new_amp, frac=r["found_frac"]))
        print(f"  {r['pair']}/{r['proto']}: amp_start {orig_amp:.4f} -> "
              f"{new_amp:.4f} nA (x{r['found_frac']})")

    # Merge into any existing manifest (keyed by pair/proto) rather than
    # overwrite -- a later --apply-file recheck must not erase the audit
    # trail of an earlier apply.
    manifest = os.path.join(SCRATCH_ROOT, "applied_patches.json")
    prior = {}
    if os.path.isfile(manifest):
        with open(manifest) as f:
            prior = {(p["pair"], p["proto"]): p for p in json.load(f)}
    for p in applied:
        prior[(p["pair"], p["proto"])] = p
    with open(manifest, "w") as f:
        json.dump(list(prior.values()), f, indent=2)
    print(f"\npatched {len(applied)} workdirs, manifest -> {manifest} "
          f"({len(prior)} total entries)")
    print("revert any one with:  mv <workdir>/prefire_simulation_config.json.orig_amp_backup "
          "<workdir>/prefire_simulation_config.json")


def _tightest_bracket(history, expected):
    """Narrowest (lo, hi) from prior history with n(lo) < expected <= n(hi).
    Assumes local monotonicity, which held for all 7 unconverged targets in
    the first run (n increased smoothly with fraction, just needed a few more
    bisection steps to land on the exact integer) -- unlike the depolarization
    -block cliffs seen elsewhere, which jump by 40+ spikes in one step and are
    NOT candidates for this retry path."""
    below = [(f, n) for f, n in history if n is not None and n < expected]
    above = [(f, n) for f, n in history if n is not None and n >= expected]
    if not below or not above:
        return FRAC_LO, FRAC_HI  # no bracket found in history, restart wide
    lo = max(below, key=lambda x: x[0])[0]
    hi = min((f for f, n in above if f > lo), default=FRAC_HI)
    return lo, hi


def _search_one_target(task, seed_lo=FRAC_LO, seed_hi=FRAC_HI, seed_history=None,
                        max_iters=MAX_ITERS):
    pair, proto, orig_n, expected = task
    trial_root = os.path.join(SCRATCH_ROOT, pair, proto)
    lo, hi = seed_lo, seed_hi
    history = list(seed_history) if seed_history else []
    found = None
    for i in range(max_iters):
        frac = round((lo + hi) / 2, 5)
        trial_dir = os.path.join(trial_root, f"amp_{frac}")
        n = _run_one_amplitude(pair, proto, frac, trial_dir)
        history.append((frac, n))
        if n is None:
            break
        if n == expected:
            found = frac
            break
        elif n > expected:
            hi = frac
        else:
            lo = frac
    return dict(pair=pair, proto=proto, orig_n=orig_n, expected=expected,
                found_frac=found, history=history,
                converged=found is not None)


def main():
    global SIM_ROOT, REFERENCE_HASH, SCRATCH_ROOT
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stdp-name", default=DEFAULT_STDP_NAME,
                    help=f"*_STDP subtree to search, e.g. L23PC_L5TTPC_STDP "
                         f"(default: {DEFAULT_STDP_NAME}). Also selects the "
                         f"scratch subdir under {SCRATCH_ROOT}, so different "
                         "cell-type sets never mix results.")
    ap.add_argument("--reference-hash", default=DEFAULT_REFERENCE_HASH,
                    help="Completed param_hash to mine guardrail_forced=True "
                         f"targets from (default: {DEFAULT_REFERENCE_HASH}, "
                         "which only exists for L5TTPC_L5TTPC_STDP -- a "
                         "different --stdp-name needs its own completed run "
                         "and its hash passed here, e.g. b8c7ff3ecf0a).")
    ap.add_argument("--pilot", type=int, default=None,
                    help="Only search the first N known-bad targets (sanity check).")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--retry-unconverged", action="store_true",
                    help="Read results.json, keep the converged entries as-is, "
                         "and re-search only the unconverged ones -- seeded "
                         "from the tightest (lo, hi) bracket already found in "
                         "their history, for --extra-iters more bisection "
                         "steps, instead of restarting from [0.5, 1.0]. Only "
                         "sound for targets that narrowed smoothly (no "
                         "depolarization-block cliff) -- check the printed "
                         "history before trusting this for a given target.")
    ap.add_argument("--extra-iters", type=int, default=6,
                    help="Bisection steps for --retry-unconverged (default 6, "
                         "resolving roughly (hi-lo)/64 of the seeded bracket).")
    ap.add_argument("--apply", action="store_true",
                    help="Write every converged found_frac into the REAL "
                         "workdir's prefire_simulation_config.json (all "
                         "converged entries in results.json must exist first "
                         "-- run the search, and --retry-unconverged if "
                         "needed, before this). Preserves the untouched "
                         "original as *.orig_amp_backup for revert. Does not "
                         "run any simulation.")
    ap.add_argument("--apply-file", default=None,
                    help="With --apply, patch from this results file instead "
                         "of results.json (e.g. a --recheck output).")
    ap.add_argument("--fit-preset", choices=sorted(FIT_PRESETS), default="cpreonly",
                    help="Which a-params/gammas drive the search (default: "
                         "cpreonly, the original validation preset). Spike "
                         "count is a-param independent far from threshold but "
                         "NOT right at a tuned threshold -- 8 of 27 amplitudes "
                         "found under cpreonly missed the target once run for "
                         "real under defit2's very different gammas. Use "
                         "--fit-preset defit2 to re-search under the "
                         "parameters you are actually about to run.")
    ap.add_argument("--only", default=None,
                    help="Comma-separated pair/protocol entries (e.g. "
                         "'185384-197412/10Hz_10ms,201812-187742/10Hz_5ms') to "
                         "re-search fresh under --fit-preset, pulling their "
                         "'expected' spike count from results.json but "
                         "ignoring its found_frac. Writes to "
                         "recheck_<preset>.json instead of results.json, so "
                         "a preset-specific fix never overwrites another "
                         "preset's validated one -- apply it separately with "
                         "--apply --apply-file <that path>.")
    args = ap.parse_args()

    global FIT_ARGS
    FIT_ARGS = FIT_PRESETS[args.fit_preset] + [
        "--tau_effca_GB_GluSynapse=278.3177658387",
        "--param_hash=calib_scan",
        f"--cpre-cpost-cache={CACHE}",
        f"--circuit-config={CIRCUIT_CONFIG}",
        "--force", "--lean",
    ]

    SIM_ROOT = sim_root_for(args.stdp_name)
    REFERENCE_HASH = args.reference_hash
    if not os.path.isdir(SIM_ROOT):
        sys.exit(f"no such results dir: {SIM_ROOT}")
    # Default stdp-name keeps the original flat path (already holds the L5
    # results); any other stdp-name gets its own subdir so results never mix.
    SCRATCH_ROOT = (os.path.join(PLASTYFIRE_ROOT, "paired_calibration")
                    if args.stdp_name == DEFAULT_STDP_NAME else
                    os.path.join(PLASTYFIRE_ROOT, "paired_calibration", args.stdp_name))

    out = os.path.join(SCRATCH_ROOT, "results.json")

    if args.only:
        base = json.load(open(out))
        wanted = set(args.only.split(","))
        todo = [r for r in base if f"{r['pair']}/{r['proto']}" in wanted]
        missing = wanted - {f"{r['pair']}/{r['proto']}" for r in todo}
        if missing:
            sys.exit(f"not found in {out}: {missing}")
        print(f"re-searching {len(todo)} target(s) fresh under --fit-preset {args.fit_preset}")
        tasks = [(r["pair"], r["proto"], r["orig_n"], r["expected"]) for r in todo]
        if args.dry_run:
            for t in tasks:
                print(f"  {t[0]}/{t[1]}  ({t[2]}/{t[3]} spikes)")
            return
        os.makedirs(SCRATCH_ROOT, exist_ok=True)
        import multiprocessing
        with multiprocessing.Pool(min(args.workers, max(1, len(tasks)))) as pool:
            results = pool.map(_search_one_target, tasks)
        for r in results:
            tag = f"amp x{r['found_frac']}" if r["converged"] else "DID NOT CONVERGE"
            print(f"{r['pair']}/{r['proto']}  {r['orig_n']}->{r['expected']} spikes  {tag}  "
                  f"steps={[(f, n) for f, n in r['history']]}")
        recheck_out = os.path.join(SCRATCH_ROOT, f"recheck_{args.fit_preset}.json")
        with open(recheck_out, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nwrote {recheck_out}")
        print(f"apply with:  python {os.path.basename(__file__)} --stdp-name {args.stdp_name} "
              f"--apply --apply-file {recheck_out}")
        return

    if args.apply:
        results = json.load(open(args.apply_file or out))
        n_bad = sum(not r["converged"] for r in results)
        if n_bad:
            sys.exit(f"{n_bad} entries in results.json are still unconverged -- "
                     f"run --retry-unconverged first, or exclude them by hand.")
        apply_patches(results)
        return

    if args.retry_unconverged:
        prev = json.load(open(out))
        done = [r for r in prev if r["converged"]]
        todo = [r for r in prev if not r["converged"]]
        print(f"kept as-is : {len(done)} already-converged entries")
        print(f"retrying   : {len(todo)} unconverged, {args.extra_iters} more steps each")
        retry_args = []
        for r in todo:
            task = (r["pair"], r["proto"], r["orig_n"], r["expected"])
            lo, hi = _tightest_bracket(r["history"], r["expected"])
            print(f"  {r['pair']}/{r['proto']}: bracket [{lo}, {hi}]")
            retry_args.append((task, lo, hi, r["history"], args.extra_iters))
        if args.dry_run:
            return
        os.makedirs(SCRATCH_ROOT, exist_ok=True)
        import multiprocessing
        with multiprocessing.Pool(min(args.workers, max(1, len(retry_args)))) as pool:
            retried = pool.starmap(_search_one_target, retry_args)
        for r in retried:
            tag = f"amp x{r['found_frac']}" if r["converged"] else "STILL DID NOT CONVERGE"
            print(f"{r['pair']}/{r['proto']}  {r['orig_n']}->{r['expected']} spikes  {tag}  "
                  f"steps={[(f, n) for f, n in r['history']]}")
        results = done + retried

    else:
        targets = find_targets()
        if args.pilot:
            targets = targets[:args.pilot]
        print(f"targets  : {len(targets)}  (from guardrail_forced=True under {REFERENCE_HASH})")
        print(f"scratch  : {SCRATCH_ROOT}")
        print(f"search   : fraction in [{FRAC_LO}, {FRAC_HI}], up to {MAX_ITERS} steps/target")
        for pair, proto, n, exp in targets[:10]:
            print(f"  {pair} / {proto}  ({n}/{exp} spikes)")
        if len(targets) > 10:
            print(f"  ... {len(targets)} total")
        if args.dry_run:
            return

        os.makedirs(SCRATCH_ROOT, exist_ok=True)
        import multiprocessing
        results = []
        with multiprocessing.Pool(min(args.workers, max(1, len(targets)))) as pool:
            for i, r in enumerate(pool.imap_unordered(_search_one_target, targets), 1):
                tag = f"amp x{r['found_frac']}" if r["converged"] else "DID NOT CONVERGE"
                print(f"[{i}/{len(targets)}] {r['pair']}/{r['proto']}  "
                      f"{r['orig_n']}->{r['expected']} spikes  {tag}  "
                      f"steps={[(f, n) for f, n in r['history']]}")
                results.append(r)

    ok = sum(r["converged"] for r in results)
    print(f"\n{ok}/{len(results)} converged to an amplitude giving exactly the target spike count")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
