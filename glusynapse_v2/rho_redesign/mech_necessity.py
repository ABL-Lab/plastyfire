"""Mechanism-necessity analysis (MECH_NECESSITY.md). Read-only; no fitter or kernel code is touched.

Part A (--floors, CSV only): groups the 54 targets by phenomenon, prints chi2 per group for every available fit tag, and the
pharmacology floors: the smallest chi2 a model can reach when it is forced to predict the same value for drug arms it cannot
tell apart (M0: control = AM251 (mglu_block) = no_block; M0 + eCB: control = no_block). Floor of a set = min over p of
sum ((m_i - p) / s_i)^2 (inverse-variance pooled p).
Part B (--pool, npz): the own spine VDCC-Ca pool of the C1 gate, c' = -c / tau_E1 + (-ica_VDCC) / i_scale (tau_E1 100 ms,
i_scale 1e-5, same recursion as scan_vgate_amp.rho_rec), read at each own pre arrival; plain, and with the t_drive 4 weight
(1 - b), b = own glutamate-bound NMDA state (tau_d_NMDA 70 ms, capped at 1). Question: can a single threshold on this pool
at own pre arrivals (eCB-LTD reusing the gate's pool) be on at the AM251-sensitive protocols and off at the AM251-insensitive
ones?
Run: sbatch, 1 CPU (cwd plastyfire, .venv). Output: logs/mech_necessity_<id>.out. MEASURED 22135544 (--floors --pool, all
POOL_SETS, 30 Zilberter pairs): 0:52, 2.09 GB MaxRSS (at the 2G limit), 85% CPU -> next full run --mem=2600M --time=0:15.
"""
import argparse, csv, glob, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE)
RES = os.path.join(HERE, "results")

# ---------------------------------------------------------------- Part A
L5, LL, ZZ = "L5->L5", "L2/3->L5", "L2/3->L2/3"
G = {  # (path, target, condition): (primary group, secondary group or "")
    (L5, "10Hz_-10ms", "control"): ("LTD", "eCB"),
    (L5, "10Hz_5ms", "control"): ("LTP", ""),
    (L5, "10Hz_10ms", "control"): ("LTP", ""),
    (L5, "10Hz_-10ms", "mglu_block"): ("eCB", "LTD"),
    (L5, "sjostrom_0.1hz_dt-10ms", "mglu_block"): ("eCB", "LTD"),
    (L5, "sjostrom_20hz_dt-10ms", "mglu_block"): ("eCB", "LTD"),
    (L5, "sjostrom_0.1hz_dt-10ms", "nmdar_block"): ("NMDA", ""),
    (L5, "sjostrom_20hz_dt-10ms", "nmdar_block"): ("NMDA", ""),
    (L5, "sjostrom_50hz_dt+10ms", "mglu_block"): ("eCB", "LTP"),
    (L5, "sjostrom_50hz_dt+10ms", "nmdar_block"): ("NMDA", ""),
    (L5, "sjostrom_0.1hz_dt-25ms", "control"): ("LTD", ""),
    (L5, "sjostrom_0.1hz_dt-120ms", "control"): ("NULL", "LTD"),
    (L5, "sjostrom_0.1hz_dt-200ms", "control"): ("NULL", "LTD"),
    (L5, "sjostrom_burst5x20hz_r50_dt-120ms", "control"): ("BURST", "LTD"),
    (L5, "sjostrom_burst5x20hz_r50_dt-200ms", "control"): ("BURST", "LTD"),
    (L5, "sjostrom_0.1hz_dt+10ms", "control"): ("LTP", ""),
    (L5, "sjostrom_0.1hz_dt-10ms", "control"): ("LTD", "eCB"),
    (L5, "sjostrom_10hz_dt+10ms", "control"): ("LTP", ""),
    (L5, "sjostrom_10hz_dt-10ms", "control"): ("LTD", ""),
    (L5, "sjostrom_20hz_dt+10ms", "control"): ("LTP", ""),
    (L5, "sjostrom_20hz_dt-10ms", "control"): ("LTD", "eCB"),
    (L5, "sjostrom_40hz_dt+10ms", "control"): ("LTP", ""),
    (L5, "sjostrom_40hz_dt-10ms", "control"): ("LTP", "BURST"),
    (L5, "sjostrom_50hz_dt+10ms", "control"): ("LTP", ""),
    (L5, "sjostrom_50hz_dt-10ms", "control"): ("LTP", "BURST"),
    (L5, "sjostrom07_step200ms_pair", "control"): ("LTP", "eCB"),
    (L5, "sjostrom07_step200ms_pair", "mglu_block"): ("eCB", "LTP"),
    (L5, "sjostrom07_step200ms_pair", "no_block"): ("NO", "LTP"),
    (L5, "sjostrom07_step200ms_pre_only", "control"): ("NULL", ""),
    (L5, "sjostrom07_step200ms_post_only", "control"): ("NULL", ""),
    (LL, "sjostrom_50hz_dt+10ms@distal", "control"): ("DIST", ""),
    (LL, "letzkus_1ap_dt+10ms", "control"): ("DIST", ""),
    (LL, "letzkus_3ap_200hz_dt+10ms@proximal", "control"): ("DIST", "LTP"),
    (LL, "letzkus_3ap_200hz_dt+10ms@distal", "control"): ("DIST", ""),
    (LL, "letzkus_3ap_200hz_dt-10ms@proximal", "control"): ("DIST", "LTD"),
    (LL, "sjostrom_50hz_dt+10ms", "control"): ("DIST", "LTP"),
    (LL, "letzkus_3ap_200hz_dt+10ms@distal", "nmdar_block"): ("NMDA", ""),
    (LL, "letzkus_nopost", "control"): ("NULL", ""),
    (LL, "letzkus_3ap_200hz_dt-500ms", "control"): ("NULL", ""),
    (ZZ, "zilberter_1ap_dt+10ms", "control"): ("LTP", ""),
    (ZZ, "zilberter_1ap_dt-10ms", "control"): ("LTD", ""),
    (ZZ, "zilberter_pre_only", "control"): ("NULL", ""),
    (ZZ, "zilberter_5ap_10hz_dt+10ms", "control"): ("LTP", "BURST"),
    (ZZ, "zilberter_5ap_20hz_dt+10ms", "control"): ("LTP", "BURST"),
    (ZZ, "zilberter_5ap_20hz_dt-10ms", "control"): ("LTD", "BURST"),
    (ZZ, "zilberter_train10_50hz_dt+4ms_last", "control"): ("BURST", ""),
    (ZZ, "zilberter_train10_50hz_dt-4ms_last", "control"): ("BURST", ""),
    (ZZ, "zilberter_train10_50hz_dt-10ms_last", "control"): ("BURST", "LTD"),
    (ZZ, "zilberter_train10_50hz_dt+5ms", "control"): ("BURST", ""),
    (ZZ, "zilberter_train10_50hz_post_only", "control"): ("NULL", ""),
    (ZZ, "zilberter_train4_50hz_dt+4ms_last", "control"): ("BURST", ""),
    (ZZ, "zilberter_train8_50hz_dt+4ms_last", "control"): ("BURST", ""),
    (ZZ, "zilberter_train10_50hz_dt+4ms_last", "mglu_block"): ("eCB", "BURST"),
    (ZZ, "zilberter_train10_50hz_dt-10ms_last", "mglu_block"): ("eCB", "BURST"),
}
GROUPS = ["LTP", "LTD", "NMDA", "eCB", "NO", "DIST", "BURST", "NULL"]
TAGS = ["C1Ajn_s5", "C1Ajz_pilot", "LM0chk", "LM0_s5", "LM0_s6", "LM1_s5", "LM1_s6",
        "LM0zchk", "LM0z_s5", "LM0z_s6", "LM1z_s5", "LM1z_s6"]
SJ04 = [("sjostrom04_dltd_pooled", "control", 0.69, 0.04), ("sjostrom04_dltd_pooled", "mglu_block", 1.06, 0.05)]


def load(tag):
    rows = {}
    for path, suf in ((L5, ""), (LL, "_l23"), (ZZ, "_l23l23")):
        f = os.path.join(RES, f"v4_{tag}{suf}.csv")
        if not os.path.isfile(f):
            continue
        for r in csv.DictReader(open(f)):
            rows[(path, r["target"], r["condition"])] = (float(r["target_mean"]), float(r["target_sem"]),
                                                         float(r["pred"]), float(r["z"]))
    return rows


def floor(ms):
    w = np.array([1 / s ** 2 for _, s in ms]); m = np.array([x for x, _ in ms])
    p = float(np.sum(w * m) / np.sum(w))
    return float(np.sum(w * (m - p) ** 2)), p


def part_a():
    data = {t: load(t) for t in TAGS}
    data = {t: d for t, d in data.items() if d}
    ref = data["C1Ajz_pilot"] if "C1Ajz_pilot" in data else data["C1Ajn_s5"]
    print("## A1. Target counts per primary group (from the C1Ajz_pilot 3-pathway score)\n")
    print("| group | L5->L5 | L2/3->L5 | L2/3->L2/3 | total |\n|---|---|---|---|---|")
    for g in GROUPS:
        c = [sum(1 for k in ref if k in G and G[k][0] == g and k[0] == p) for p in (L5, LL, ZZ)]
        print(f"| {g} | {c[0]} | {c[1]} | {c[2]} | {sum(c)} |")
    miss = [k for k in ref if k not in G]
    print(f"\nunmapped rows: {miss}\n")

    print("## A2. Pharmacology floors (data only, no model)\n")
    print("| path | target | arms forced equal | data (mean +- sem) | pooled p | floor chi2 |\n|---|---|---|---|---|---|")
    tgts = sorted({(k[0], k[1]) for k in ref})
    tot = {"M0": 0.0, "M0+eCB": 0.0}
    for path, tg in tgts:
        arms = {c: ref[(path, tg, c)][:2] for c in ("control", "mglu_block", "no_block") if (path, tg, c) in ref}
        if len(arms) < 2:
            continue
        for lvl, keep in (("M0", ("control", "mglu_block", "no_block")), ("M0+eCB", ("control", "no_block"))):
            sel = [arms[c] for c in keep if c in arms]
            if len(sel) < 2:
                continue
            f, p = floor(sel)
            tot[lvl] += f
            ds = ", ".join(f"{c} {arms[c][0]:.3f}+-{arms[c][1]:.3f}" for c in keep if c in arms)
            print(f"| {path} | {tg} | {lvl}: {'='.join(c for c in keep if c in arms)} | {ds} | {p:.3f} | {f:.2f} |")
    f, p = floor([(m, s) for _, _, m, s in SJ04])
    print(f"| {L5} | sjostrom04 dLTD (planned, pooled) | M0: control=mglu_block | 0.69+-0.04, 1.06+-0.05 | {p:.3f} | {f:.2f} |")
    print(f"\nsum of floors over the 54 scored targets: M0 {tot['M0']:.2f}, M0+eCB {tot['M0+eCB']:.2f}; sj04 adds {f:.2f} "
          f"to M0 and M0+eCB-without-sj04-drive.\n")
    nm = [(k, v) for k, v in ref.items() if k[2] == "nmdar_block"]
    print("nmdar_block rows (rho, dpre frozen in every model, so constant): " +
          ", ".join(f"{k[1]} pred {v[2]:.3f} z {v[3]:+.2f}" for k, v in nm) +
          f"; constant chi2 {sum(v[3] ** 2 for _, v in nm):.2f}\n")

    print("## A3. chi2 per primary group, per fit tag (n = targets scored in that tag)\n")
    print("| tag | " + " | ".join(GROUPS) + " | total (n) |\n|---|" + "---|" * (len(GROUPS) + 1))
    for t, d in data.items():
        cells = []
        for g in GROUPS:
            zs = [v[3] for k, v in d.items() if k in G and G[k][0] == g]
            cells.append(f"{sum(z * z for z in zs):.1f} ({len(zs)})" if zs else "-")
        print(f"| {t} | " + " | ".join(cells) + f" | {sum(v[3] ** 2 for v in d.values()):.1f} ({len(d)}) |")
    print("\n## A4. Per-target predictions (pred / z) for every tag\n")
    hdr = list(data)
    print("| path | target | cond | group | data | " + " | ".join(hdr) + " |\n|---|---|---|---|---|" + "---|" * len(hdr))
    for k in sorted(G, key=lambda k: (GROUPS.index(G[k][0]), k)):
        if k not in ref:
            continue
        cells = [f"{data[t][k][2]:.2f} / {data[t][k][3]:+.1f}" if k in data[t] else "-" for t in hdr]
        print(f"| {k[0]} | {k[1]} | {k[2]} | {G[k][0]}{('/' + G[k][1]) if G[k][1] else ''} | "
              f"{ref[k][0]:.2f}+-{ref[k][1]:.2f} | " + " | ".join(cells) + " |")
    print()


# ---------------------------------------------------------------- Part B
POOL_SETS = [  # (dir, protocol, basis, class)
    ("ebner_delta-prefire-vca", "sjostrom_0.1hz_dt-10ms", "L5", "eCB+ (AM251 abolishes)"),
    ("ebner_delta-prefire-vca", "sjostrom_20hz_dt-10ms", "L5", "eCB+ (AM251 abolishes)"),
    ("markram_delta-prefire-vca", "10Hz_-10ms", "L5", "eCB+ (AM251 abolishes)"),
    ("sj07_delta-prefire-vca", "sjostrom07_step200ms_pair", "L5", "eCB+ (AM251 enhances LTP)"),
    ("sj04_delta-prefire-vca", "sjostrom04_*", "L5", "eCB+ (AM251 abolishes dLTD)"),
    ("ebner_delta-prefire-vca", "sjostrom_50hz_dt+10ms", "L5", "eCB- (AM251 no effect)"),
    ("zilberter_l23l23_delta-prefire-vseg-rs", "zilberter_train10_50hz_dt-10ms_last", "ZZ", "eCB- (AM251 no effect)"),
    ("zilberter_l23l23_delta-prefire-vseg-rs", "zilberter_train10_50hz_dt+4ms_last", "ZZ", "eCB n.s. (AM251 1.73 vs 1.49)"),
    ("sj03_delta-prefire-vca", "sjostrom_0.1hz_dt-25ms", "L5", "LTD, drug untested"),
    ("ebner_delta-prefire-vca", "sjostrom_10hz_dt-10ms", "L5", "LTD, drug untested"),
    ("sj03r50_delta-prefire-vca", "sjostrom_burst5x20hz_r50_dt-120ms", "L5", "LTD, drug untested"),
    ("sj03r50_delta-prefire-vca", "sjostrom_burst5x20hz_r50_dt-200ms", "L5", "LTD, drug untested"),
    ("zilberter_l23l23_delta-prefire-vseg-rs", "zilberter_1ap_dt-10ms", "ZZ", "LTD, drug untested"),
    ("zilberter_l23l23_delta-prefire-vseg-rs", "zilberter_1ap_dt+10ms", "ZZ", "LTD (+10), drug untested"),
    ("sj03_delta-prefire-vca", "sjostrom_0.1hz_dt-120ms", "L5", "no change"),
    ("sj03_delta-prefire-vca", "sjostrom_0.1hz_dt-200ms", "L5", "no change"),
    ("ebner_delta-prefire-vca", "sjostrom_0.1hz_dt+10ms", "L5", "no change"),
    ("ebner_delta-prefire-vca", "sjostrom_10hz_dt+10ms", "L5", "LTP"),
    ("ebner_delta-prefire-vca", "sjostrom_50hz_dt-10ms", "L5", "LTP"),
]
THS = [0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]   # in units of theta_VDCC (C1Ajn_s5)


def pool_at_arrivals(t, vd, ta, ka, tau, isc, tau_d):
    """vd (T,) own -ica_VDCC; ta arrival times, ka arrival grid indices. Returns (plain, weighted) pool at each arrival.
    c[k] = sum_{j<k} bV[j] vd[j] exp(-(t[k] - t[j+1]) / tau), bV = tau (1 - exp(-h / tau)) / isc (= kernel recursion)."""
    h = np.diff(t); bV = tau * (1.0 - np.exp(-h / tau)) / isc
    out = np.zeros((len(ta), 2))
    for n, (a, k) in enumerate(zip(ta, ka)):
        j0 = int(np.searchsorted(t, a - 1200.0))
        if k <= j0:
            continue
        j = np.arange(j0, k)
        dec = np.exp(-(t[k] - t[j + 1]) / tau)
        x = bV[j] * vd[j] * dec
        tj = t[j]; b = np.zeros(len(j))
        for a2 in ta[(ta > a - 1500.0) & (ta <= t[k])]:
            m = tj >= a2 - 1e-9
            b[m] += np.exp(-(tj[m] - a2) / tau_d)
        out[n] = (x.sum(), (x * (1.0 - np.minimum(b, 1.0))).sum())
    return out


def part_b(fit, max_pairs_zz, only=""):
    sys.path.insert(0, V2)
    import batch_v2
    import model_v2 as MV
    from batch_v2 import BatchV2
    fj = json.load(open(fit))
    P = {**json.loads(fj["args"]["filters"]), **json.loads(fj["args"].get("set") or "{}"), **fj["pre"]}
    thV = float(P["theta_V"]); tau = float(P["tau_E1"]); isc = float(P["i_scale"]); tau_d = MV.DEFAULTS["tau_d_NMDA"]
    print(f"## B. Own spine VDCC-Ca pool at own pre arrivals (fit {os.path.basename(fit)}: theta_VDCC {thV:.3f}, "
          f"tau_E1 {tau}, i_scale {isc})\n")
    print("Per synapse: pool at each arrival, in units of theta_VDCC; 'max' = max over arrivals, 'frac>th' = fraction of\n"
          "arrivals above th, averaged over synapses. w = (1 - b) weighted (VDCC Ca before own glutamate).\n")
    l5pairs = set(open(os.path.join(V2, "subset24_pairs.txt")).read().strip().split(","))
    bases = {"L5": batch_v2.BASIS_DIR, "ZZ": os.path.join(batch_v2.ROOT, "basis_results_edges_zilberter_l23l23_delta_rs")}
    rows = []
    for dname, proto, bk, cls in [s for s in POOL_SETS if only in s[1] or only in s[3]]:
        batch_v2.BASIS_DIR = bases[bk]
        fs = sorted(glob.glob(os.path.join(V2, "extracted", dname, f"*__{proto}.npz")))
        if bk == "L5":
            fs = [f for f in fs if os.path.basename(f).split("__")[0] in l5pairs]
        else:
            fs = fs[:max_pairs_zz]
        by = {}
        for f in fs:
            pair, pr = os.path.basename(f)[:-4].split("__")
            B = BatchV2([os.path.dirname(f)], protocols=[pr], pairs={pair}, fast=False, signals=("vdcc",), verbose=False)
            if not B.recs or B.recs[0]["vdcc"] is None:
                continue
            r = B.recs[0]; t = np.asarray(r["t"], np.float64)
            delay = MV.edge_params(r["syn"])["delay"]
            for i in range(len(r["syn"])):
                ta = np.asarray(r["prespikes"], np.float64) + delay[i]
                if len(ta) == 0:
                    continue
                pv = pool_at_arrivals(t, np.asarray(r["vdcc"][i], np.float64), ta, np.asarray(r["arr"][i]), tau, isc,
                                      tau_d) / thV
                by.setdefault(pr, []).append(pv)
            del B
        for pr, L in by.items():
            mx = np.array([p.max(axis=0) for p in L])                       # (n_syn, 2)
            fr = {th: np.mean([np.mean(p[:, 1] > th) for p in L]) for th in THS}
            frp = {th: np.mean([np.mean(p[:, 0] > th) for p in L]) for th in THS}
            q = lambda x: "/".join(f"{v:.3g}" for v in np.percentile(x, [25, 50, 75]))
            rows.append((pr, cls, len(L), q(mx[:, 0]), q(mx[:, 1]), fr, frp))
            print(f"{pr}: {len(L)} syn done", flush=True)
    print("\n| protocol | class | n syn | max pool plain q25/50/75 | max pool w q25/50/75 | " +
          " | ".join(f"w frac>{th}" for th in THS) + " |\n|---|---|---|---|---|" + "---|" * len(THS))
    for pr, cls, n, a, b, fr, _ in rows:
        print(f"| {pr} | {cls} | {n} | {a} | {b} | " + " | ".join(f"{fr[th]:.2f}" for th in THS) + " |")
    print("\n| protocol | class | " + " | ".join(f"plain frac>{th}" for th in THS) + " |\n|---|---|" + "---|" * len(THS))
    for pr, cls, n, a, b, _, frp in rows:
        print(f"| {pr} | {cls} | " + " | ".join(f"{frp[th]:.2f}" for th in THS) + " |")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--floors", action="store_true")
    ap.add_argument("--pool", action="store_true")
    ap.add_argument("--fit", default=os.path.join(RES, "v4_C1Ajn_s5.json"))
    ap.add_argument("--max-pairs-zz", type=int, default=30)
    ap.add_argument("--only", default="", help="keep POOL_SETS whose protocol or class contains this")
    a = ap.parse_args()
    if a.floors:
        part_a()
    if a.pool:
        part_b(a.fit, a.max_pairs_zz, a.only)
