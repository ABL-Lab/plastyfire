"""
STDP at +/-10 ms: b78bbff emodel vs production, on the SHARED pair set.

b78bbff (2025-10-17) is the emodel behind the Jun 30 figure: apical gNaTgbar_NaTg
0.04, basal 0.025, geom_nseg_fixed(20). Production (912e82a) halved apical Na to
0.019856 at ec0d47c and moved to nseg 40. Everything else -- gen-4 fit params,
edges.h5, protocols -- is identical between arms, so the emodel is the only variable.

Restricted to pairs present in BOTH caches (51 of them): b78bbff only produced a
usable single-AP c_post on 52/100 pairs (those cells fire doublets at the higher Na
density), so comparing full arms would confound the emodel change with a different
pair sample.

Each arm is scored on its OWN basis -- basis_b78bbff/ vs
basis_results_edges_ion_channels/ -- since basis EPSPs are measured on the emodel.
"""

import os
import pickle
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plastyfire.evaluator_edges import compute_epsp_from_basis, _load_basis  # noqa: E402

RESULTS_DIR = os.path.join(
    HERE, "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations")
HASH = "0ce64fa83b85"
PROTO_TO_DT = {"10Hz_10ms": 10, "10Hz_-10ms": -10}
INVITRO = {-10: (0.7922, 0.0259), 10: (1.2013, 0.0626)}

ARMS = {
    "b78bbff": dict(tag="b78bbff_", basis="basis_b78bbff",
                    cache="cpre_cpost_cache/b78bbff_tau278.pkl",
                    label="b78bbff (Na 0.04, nseg 20)", color="#1565c0"),
    "production": dict(tag="", basis="basis_results_edges_ion_channels",
                       cache="cpre_cpost_cache/ion_channels_tau278.pkl",
                       label="production (Na 0.0199, nseg 40)", color="#2e7d32"),
}


def ratio(arm, pair, proto):
    a = ARMS[arm]
    pkl = os.path.join(RESULTS_DIR, pair, proto,
                       "simulation_edges_%s%s.pkl" % (a["tag"], HASH))
    if not os.path.exists(pkl):
        return None
    try:
        raw = pickle.load(open(pkl, "rb"))
        pre, post = pair.split("-")
        bdf = _load_basis(os.path.join(HERE, a["basis"]), pre, post)
        b_m, b_s = compute_epsp_from_basis(bdf, list(raw["initial_rho"]))
        a_m, _ = compute_epsp_from_basis(bdf, list(raw["final_rho"]))
        if b_m == 0:
            return None
        return (a_m / b_m) * (1.0 + min((b_s / b_m) ** 2, 0.25))
    except Exception:
        return None


def main():
    caches = {k: set(pickle.load(open(os.path.join(HERE, v["cache"]), "rb")).keys())
              for k, v in ARMS.items()}
    shared = caches["b78bbff"] & caches["production"]
    print("b78bbff cache : %d pairs" % len(caches["b78bbff"]))
    print("production    : %d pairs" % len(caches["production"]))
    print("shared        : %d pairs\n" % len(shared))
    shared_names = {"%d-%d" % p for p in shared}

    pairs = sorted(p for p in os.listdir(RESULTS_DIR)
                   if re.fullmatch(r"\d+-\d+", p) and p in shared_names)

    data = {k: {d: [] for d in PROTO_TO_DT.values()} for k in ARMS}
    paired = {d: [] for d in PROTO_TO_DT.values()}   # (b78, prod) on same pair
    for pair in pairs:
        for proto, dt in PROTO_TO_DT.items():
            rb = ratio("b78bbff", pair, proto)
            rp = ratio("production", pair, proto)
            if rb is not None:
                data["b78bbff"][dt].append(rb)
            if rp is not None:
                data["production"][dt].append(rp)
            if rb is not None and rp is not None:
                paired[dt].append((rb, rp))

    dts = sorted(PROTO_TO_DT.values())
    print("Mean EPSP ratio (+/- SEM):")
    for k in ARMS:
        parts = []
        for d in dts:
            v = np.array(data[k][d])
            parts.append("dt=%+3d: %.4f +/- %.4f (n=%d)" % (
                d, v.mean(), v.std() / np.sqrt(len(v)), len(v)) if len(v) else "--")
        print("  %-30s %s" % (ARMS[k]["label"], "   ".join(parts)))
    print("  %-30s %s" % ("in vitro (Markram 1997)",
          "   ".join("dt=%+3d: %.4f +/- %.4f" % (d, *INVITRO[d]) for d in dts)))

    print("\nPaired (same pair, both arms):")
    for d in dts:
        if paired[d]:
            arr = np.array(paired[d])
            diff = arr[:, 0] - arr[:, 1]
            print("  dt=%+3d  n=%d  b78=%.4f  prod=%.4f  mean diff=%+.4f +/- %.4f" % (
                d, len(arr), arr[:, 0].mean(), arr[:, 1].mean(),
                diff.mean(), diff.std() / np.sqrt(len(diff))))

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for k in ARMS:
        m = [np.mean(data[k][d]) for d in dts]
        s = [np.std(data[k][d]) / np.sqrt(len(data[k][d])) for d in dts]
        ax.errorbar(dts, m, yerr=s, marker="o", capsize=4, lw=2.2,
                    color=ARMS[k]["color"], label=ARMS[k]["label"])
    ax.errorbar(dts, [INVITRO[d][0] for d in dts], yerr=[INVITRO[d][1] for d in dts],
                marker="^", ls="--", capsize=4, lw=2, color="#ef6c00",
                label="in vitro (Markram 1997)")
    ax.axhline(1.0, color="grey", ls="--", lw=0.8)
    ax.axvline(0, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("dt (ms)")
    ax.set_ylabel("Mean EPSP ratio")
    ax.set_title("Old emodel (b78bbff) vs production\n%d shared pairs, gen-4 params" %
                 len(pairs))
    ax.set_xticks(dts)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(HERE, "stdp_b78bbff_vs_prod.png")
    plt.savefig(out, dpi=200)
    print("\nSaved -> %s" % out)


if __name__ == "__main__":
    main()
