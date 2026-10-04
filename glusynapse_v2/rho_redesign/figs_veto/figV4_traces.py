"""Inputs of figV4_why_ecb.png: why the eCB step exists (0.1 Hz single pairings) and why it then needs the veto (10 Hz).

Protocols (L5 -> L5): Sjostrom 0.1 Hz -10 / +10, Markram 10 Hz -10 / +10.
Fits: B4 (u97, k_V 0.75), KOE = B4 with A_eCB 0 (eCB step off), KOV = B4 with veto_peak_k 1e30 (veto off),
CHR = Chindemi refit (CHR_s5). No refit; nothing in the rule changes.
1. Synapse: the figV3 synapse (rho0 = 0, pair 181455-195199), traced under B4: cai, pool, events.
2. Population: every synapse of every pair, control-lane d and rho on a common grid; per-pair EPSP ratio (ratio_hook).
-> /scratch/dhuruva/figs_veto/figV4.npz, figV4_meta.json, figV4_outcomes.csv

  python glusynapse_v2/rho_redesign/figs_veto/figV4_traces.py
(env as figV_traces.py). Size: figV3_traces 22480554 (2 protocols x 3 fits, 86 s); figV_traces 22479950 2:21, 3.60 GB
-> 4 x 4 here: 4500M, 0:15. MEASURED 22480956 (traces + plot) 2:20, 4.39 GB -> next run 5500M, 0:15.
"""
import copy, gc, json, os, sys, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figV_traces as FV                                  # noqa: E402  (imports tv12_lat as FV.T)

T = FV.T
WORK = "/scratch/dhuruva/figs_veto"
PROTOS = dict(s01m="sjostrom_0.1hz_dt-10ms", s01p="sjostrom_0.1hz_dt+10ms", m10="10Hz_-10ms", p10="10Hz_10ms")
FITS = ("B4", "KOE", "KOV", "CHR")
PATH = "L5"
NG = 1500


def koe(fit):
    f = copy.deepcopy(fit)
    f["pre"]["A_eCB"] = 0.0
    fl = json.loads(f["args"]["filters"]); fl["A_eCB"] = 0.0; f["args"]["filters"] = json.dumps(fl)
    f["variant"] = dict(name="KOE", base=FV.B4J, filter_changes={"A_eCB": 0.0}, refit=False)
    return f


def main():
    import argparse
    T.A = argparse.Namespace(out=WORK); T.SEL, T.STATS, T.CHOICE, T.VFILES = [], {}, {}, {}
    b4 = json.load(open(FV.B4J))
    T._FJ.update(B4=b4, KOE=koe(b4), KOV=FV.kov(b4), CHR=json.load(open(FV.CHRJ)))
    t0 = time.time()
    ch = json.load(open(os.path.join(WORK, "figV3_meta.json")))["choice"]
    pr, syn = ch["pair"], int(ch["syn"])
    val = pd.read_csv(FV.GPUVAL["B4"][0]); val = val[val.condition == "control"].set_index("target")
    pp = pd.read_csv(FV.GPUVAL["B4"][1])
    pp = pp[(pp.pathway == "l5") & (pp.condition == "control")]

    pkt, meta = {}, dict(choice=dict(pair=pr, syn=syn), protos=PROTOS, ev_fields=list(T.EV_FIELDS), traces={},
                         port_fail=[])
    out = []
    for k, q in PROTOS.items():
        pairs = sorted(set(pp[pp.target == q].pair))
        D = {fk: [] for fk in FITS}; R = {fk: [] for fk in FITS}; RAT = {fk: {} for fk in FITS}; R0 = []
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
                RAT[fk][p] = float(np.nanmean(T.ratio_hook(r, args, info, ref)["ratio_ctl"]))
                TR, EV, NEV, tr = T.trace_run(args, n, Tn, list(range(n)), info["cnt"])
                if not T.check(f"{k} {p} {fk}", ref, tr, TR, np.arange(n))["ok"]:
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
                    pkt[f"{k}__pool"] = TR[i, T.F_["pool"]].astype(np.float32)
                    pkt[f"{k}__EV"] = EV[i, :int(NEV[i])]
                    meta["traces"][k] = dict(i=i, uE=float(info["uE"][i]), thE=float(info["thE"][i]),
                                             vt0=info["vt0"], vtv=info["vtv"], rho0=float(r["rho0"][i]))
                del TR, EV
            R0 += [float(v) for v in r["rho0"]]
            T._REC.pop((PATH, p, q), None); del r; gc.collect()
            T.drop_recs(); gc.collect()
            print(f"  {k} {p}: done ({time.time() - t0:.0f} s)", flush=True)
        pkt[f"{k}__grid"] = grid; pkt[f"{k}__rho0"] = np.array(R0)
        for fk in FITS:
            pkt[f"{k}__{fk}__D"] = np.array(D[fk], np.float32); pkt[f"{k}__{fk}__R"] = np.array(R[fk], np.float32)
            x = pd.Series(RAT[fk]).dropna()
            out.append(dict(proto=k, target=q, fit=fk, mean=x.mean(), sem=x.std(ddof=1) / np.sqrt(len(x)), n=len(x),
                            data_mean=val.loc[q, "target_mean"], data_sem=val.loc[q, "target_sem"],
                            gpu_pred=val.loc[q, "pred"] if fk == "B4" else np.nan))
            print(f"POP {k} {fk}: ratio {x.mean():.3f} +/- {out[-1]['sem']:.3f} (n {len(x)}); d end "
                  f"{np.array(D[fk])[:, -1].mean():+.3f}; rho end {np.array(R[fk])[:, -1].mean():.3f}", flush=True)
        if k in meta["traces"]:
            E = pkt[f"{k}__EV"]; ix = {f: T.EV_FIELDS.index(f) for f in T.EV_FIELDS}; s = meta["traces"][k]
            for j in range(min(4, len(E))):
                print(f"  {k} arrival {j + 1}: t {E[j, ix['t_a']] + pkt[f'{k}__t'][0] - pkt[f'{k}__pre'][0]:8.1f} ms "
                      f"pool/P1 {E[j, ix['trigval']] / s['uE']:.3f} trig {int(E[j, ix['trig']])} "
                      f"veto {int(E[j, ix['veto']])} step {int(E[j, ix['step']])}", flush=True)
    oc = pd.DataFrame(out)
    oc.to_csv(os.path.join(WORK, "figV4_outcomes.csv"), index=False)
    print(oc.to_string(index=False, float_format="%.3f"), flush=True)
    np.savez_compressed(os.path.join(WORK, "figV4.npz"), **pkt)
    json.dump(meta, open(os.path.join(WORK, "figV4_meta.json"), "w"), indent=1, default=str)
    print(f"\nwrote {WORK}/figV4.npz ({time.time() - t0:.0f} s); port failures {meta['port_fail'] or 'none'}", flush=True)
    if meta["port_fail"]:
        sys.exit(3)


if __name__ == "__main__":
    main()
