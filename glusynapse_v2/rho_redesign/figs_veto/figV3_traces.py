"""Inputs of figV3_why_veto.png: Markram 10 Hz -10 vs +10, one rho0 = 0 synapse + population time courses.

Same machinery and parameter sets as figV_traces.py (B4 = u97 k0.75, KOV = B4 with veto_peak_k 1e30, CHR = CHR_s5);
nothing in the rule or the parameters is changed.
1. Synapse (rho0 = 0): among rho0 = 0 synapses present in both records that step at 10 Hz -10 and have a vetoed
   trigger at 10 Hz +10 under B4 (figV_pop.csv), the one closest to the class median of rho_end under B4 at both
   protocols (sum of |rho_end - median| / IQR). The figV_traces contrast |d_KOV - d_B4| is 0.29 for every candidate,
   so the old median rule cannot choose. Traced under B4: t, cai (recovered from effcai), c*, pool, EV.
2. Population: every synapse of every pair of the two targets, traced under B4 / KOV / CHR; control-lane d and rho
   resampled onto a common grid (s from the first pre spike, NG points); saved per synapse with rho0 and pair.
   Port check (traced end state == v12_rec) on every record.
-> /scratch/dhuruva/figs_veto/figV3.npz, figV3_meta.json

  python glusynapse_v2/rho_redesign/figs_veto/figV3_traces.py
(env as figV_traces.py). Size: MEASURED figV_traces 22479950 2:21, 3.60 GB -> 4500M, 0:15.
"""
import gc, json, os, sys, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figV_traces as FV                                  # noqa: E402  (imports tv12_lat as FV.T)

T = FV.T
WORK = "/scratch/dhuruva/figs_veto"
PROTOS = dict(m10="10Hz_-10ms", p10="10Hz_10ms")
FITS = ("B4", "KOV", "CHR")
PATH = "L5"
NG = 1500


def main():
    import argparse
    T.A = argparse.Namespace(out=WORK); T.SEL, T.STATS, T.CHOICE, T.VFILES = [], {}, {}, {}
    b4 = json.load(open(FV.B4J))
    T._FJ.update(B4=b4, KOV=FV.kov(b4), CHR=json.load(open(FV.CHRJ)))
    t0 = time.time()
    df = pd.read_csv(os.path.join(WORK, "figV_pop.csv"))

    # ---- verify the population numbers quoted by the coordinator
    print("figV_pop.csv end states (mean over synapses) by rho0 class:")
    for fk in FITS:
        for k in PROTOS:
            x = df[(df.proto == k) & (df.fit == fk)]
            for r0 in (0.0, 1.0):
                y = x[x.rho0 == r0]
                print(f"  {fk:3s} {k}: rho0 {r0:.0f} n {len(y):3d}  rho_end {y.rho.mean():.3f}  d_eCB {y.d.mean():.3f}  "
                      f"control d {y.dctl.mean():.3f}  steps/syn {y['ne'].mean():.1f}  vetoed/syn {y['nv'].mean():.1f}", flush=True)
            print(f"  {fk:3s} {k}: all n {len(x)}  rho_end {x.rho.mean():.3f}  control d {x.dctl.mean():.3f}", flush=True)

    # ---- 1. synapse choice (rho0 = 0)
    def sel(k, fk, tag):
        x = df[(df.proto == k) & (df.fit == fk)].set_index(["pair", "syn"])
        return x.drop(columns=["proto", "fit"]).add_suffix("_" + tag)
    w = sel("m10", "B4", "m").join(sel("p10", "B4", "p"), how="inner").reset_index()
    c = w[(w.rho0_m == 0) & (w.ne_m > 0) & (w.nv_p > 0)].copy()
    print(f"\nrho0 = 0 candidates (step at -10, vetoed trigger at +10): {len(c)} of {int((w.rho0_m == 0).sum())} rho0 = 0 "
          f"synapses in both records", flush=True)
    if c.empty:
        c = w[w.rho0_m == 0].copy(); print("  relaxed to all rho0 = 0 synapses", flush=True)
    sc = 0.0
    for col in ("rho_m", "rho_p"):
        iqr = max(float(c[col].quantile(.75) - c[col].quantile(.25)), 1e-3)
        sc = sc + (c[col] - c[col].median()).abs() / iqr
    c["score"] = sc
    row = c.sort_values("score").iloc[0]
    pr, syn = row["pair"], int(row["syn"])
    choice = dict(pair=pr, syn=syn, n_cand=int(len(c)), rule="rho0 = 0, steps at 10Hz -10 and vetoed trigger at 10Hz +10 "
                  "under B4; closest to the class median rho_end at both protocols (sum |rho_end - median| / IQR)",
                  med_rho_m=float(c.rho_m.median()), med_rho_p=float(c.rho_p.median()),
                  rho_m=float(row.rho_m), rho_p=float(row.rho_p))
    print(f"CHOSEN {json.dumps(choice)}", flush=True)

    pkt, meta = {}, dict(choice=choice, protos=PROTOS, ev_fields=list(T.EV_FIELDS), traces={}, port_fail=[])
    # ---- 2. population + the chosen synapse, one record at a time
    pp = pd.read_csv(FV.GPUVAL["B4"][1])
    pp = pp[(pp.pathway == "l5") & (pp.condition == "control")]
    for k, q in PROTOS.items():
        pairs = sorted(set(pp[pp.target == q].pair))
        D = {fk: [] for fk in FITS}; R = {fk: [] for fk in FITS}; R0, PRS = [], []
        grid = None
        for p in pairs:
            if not T.has_rec(PATH, p, q):
                print(f"  {k} {p}: no record", flush=True); continue
            spv, *_ = T.pair_spk(T.FJ("B4"), PATH, p, [q])
            r = T.load_rec(PATH, p, q)
            n, Tn = r["effcai"].shape
            t = np.asarray(r["t"], np.float64)
            pre = np.sort(np.asarray(r["prespikes"], float).ravel())
            xs = (t - pre[0]) / 1000.0
            if grid is None:
                grid = np.linspace(0.0, xs[-1], NG)
            for fk in FITS:
                args, info = T.prep(T.FJ(fk), PATH, r, spv)
                ref = T.ref_run(args, n)
                TR, EV, NEV, tr = T.trace_run(args, n, Tn, list(range(n)), info["cnt"])
                chk = T.check(f"{k} {p} {fk}", ref, tr, TR, np.arange(n))
                if not chk["ok"]:
                    meta["port_fail"].append(f"{k} {p} {fk}")
                dctl, _ = T.lanes(info, TR[:, T.F_["d_eCB"]], TR[:, T.F_["d_NO"]])
                for i in range(n):
                    D[fk].append(np.interp(grid, xs, dctl[i]))
                    R[fk].append(np.interp(grid, xs, TR[i, T.F_["rho"]]))
                if p == pr and fk == "B4":
                    i = int(np.flatnonzero(np.asarray(r["syn"]) == syn)[0])
                    pkt[f"{k}__t"] = t; pkt[f"{k}__pre"] = pre
                    pkt[f"{k}__post"] = np.sort(np.asarray(r["postspikes"], float).ravel())
                    pkt[f"{k}__cai"] = FV.cai_from_effcai(r["effcai"][i], t).astype(np.float32)
                    for f in ("c", "pool", "vdcc"):
                        pkt[f"{k}__{f}"] = TR[i, T.F_[f]].astype(np.float32)
                    pkt[f"{k}__EV"] = EV[i, :int(NEV[i])]
                    meta["traces"][k] = dict(i=i, uE=float(info["uE"][i]), thE=float(info["thE"][i]),
                                             thP=float(info["thP"][i]), td=float(info["td"][i]), tp=float(info["tp"][i]),
                                             vt0=info["vt0"], vtv=info["vtv"], rho0=float(r["rho0"][i]),
                                             rho_end=float(ref["rho"][i]), d_end=float(ref["d"][i]))
                del TR, EV
            R0 += [float(v) for v in r["rho0"]]; PRS += [p] * n
            T._REC.pop((PATH, p, q), None); del r; gc.collect()
            T.drop_recs(); gc.collect()
            print(f"  {k} {p}: done ({time.time() - t0:.0f} s)", flush=True)
        pkt[f"{k}__grid"] = grid; pkt[f"{k}__rho0"] = np.array(R0)
        for fk in FITS:
            pkt[f"{k}__{fk}__D"] = np.array(D[fk], np.float32); pkt[f"{k}__{fk}__R"] = np.array(R[fk], np.float32)
        r0 = np.array(R0)
        for fk in FITS:
            Dm, Rm = np.array(D[fk]), np.array(R[fk])
            print(f"POP {k} {fk}: n {len(r0)}; control d end mean {Dm[:, -1].mean():.3f}; rho end rho0=0 "
                  f"{Rm[r0 == 0, -1].mean():.3f} (n {int((r0 == 0).sum())}), rho0=1 {Rm[r0 == 1, -1].mean():.3f} "
                  f"(n {int((r0 == 1).sum())})", flush=True)
        if k in meta["traces"]:
            E = pkt[f"{k}__EV"]; ix = {f: T.EV_FIELDS.index(f) for f in T.EV_FIELDS}; s = meta["traces"][k]
            cai, t = pkt[f"{k}__cai"], pkt[f"{k}__t"]
            print(f"\nchosen {k}: rho {s['rho0']:.0f} -> {s['rho_end']:.3f}, d_eCB end {s['d_end']:.3f}", flush=True)
            for j in range(min(5, len(E))):
                a = E[j, ix["t_a"]] + t[0]
                m = (t >= a - 2) & (t < a + 60)
                p0 = pkt[f"{k}__pre"][0]
                print(f"  arrival {j + 1}: t {a - p0:7.1f} ms  pool/P1 "
                      f"{E[j, ix['trigval']] / s['uE']:.3f}  trig {int(E[j, ix['trig']])} veto {int(E[j, ix['veto']])} "
                      f"step {int(E[j, ix['step']])}  cai peak in [t_a-2, t_a+60) {np.nanmax(cai[m]):.3f} uM", flush=True)
    np.savez_compressed(os.path.join(WORK, "figV3.npz"), **pkt)
    json.dump(meta, open(os.path.join(WORK, "figV3_meta.json"), "w"), indent=1, default=str)
    print(f"\nwrote {WORK}/figV3.npz ({time.time() - t0:.0f} s); port failures {meta['port_fail'] or 'none'}", flush=True)
    if meta["port_fail"]:
        sys.exit(3)


if __name__ == "__main__":
    main()
