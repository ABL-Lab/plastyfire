"""Inputs of figV5_chindemi_amount.png: why the Chindemi refit cannot fit Sjostrom 2001 -10 / +10 at 0.1, 10, 20 Hz.

Fit: CHR = Chindemi refit (CHR_s5), L5 -> L5, control. Nothing changed.
Per synapse: time c* spends in the depression zone [theta_d, theta_p) and above theta_p over the whole induction;
per pair: EPSP ratio (ratio_hook). One synapse (the figV3/V4 synapse if it has the record): c* trace, theta_d, theta_p.
-> /scratch/dhuruva/figs_veto/figV5.npz, figV5_meta.json, figV5_zones.csv, figV5_outcomes.csv

  python glusynapse_v2/rho_redesign/figs_veto/figV5_traces.py
(env as figV_traces.py). Size: figV4 22480956 (16 runs incl. 0.1 Hz) 2:20, 4.39 GB -> 5500M, 0:15.
MEASURED 22482004 (traces + plot) 3:58, 5.36 GB -> next run 6700M, 0:15.
MEASURED 22487045 (FREQS 5, traces only) 6:02, 6.53 GB -> next 5-freq run 8200M, 0:15.
"""
import gc, json, os, sys, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figV_traces as FV                                  # noqa: E402

T = FV.T
WORK = "/scratch/dhuruva/figs_veto"
FREQS = os.environ.get("FREQS", "0.1,10,20").split(",")
TAG = os.environ.get("FIGTAG", "figV5")
PROTOS = {f"{f}{s}": f"sjostrom_{f}hz_dt{s}10ms" for f in FREQS for s in ("-", "+")}
PATH = "L5"
EBNER = "/lustre09/project/6070394/dhuruva/plastyfire/ebner/ebner_targets.csv"


def main():
    import argparse
    T.A = argparse.Namespace(out=WORK); T.SEL, T.STATS, T.CHOICE, T.VFILES = [], {}, {}, {}
    T._FJ.update(CHR=json.load(open(FV.CHRJ)))
    t0 = time.time()
    ch = json.load(open(os.path.join(WORK, "figV3_meta.json")))["choice"]
    pp = pd.read_csv(FV.GPUVAL["B4"][1])
    pp = pp[(pp.pathway == "l5") & (pp.condition == "control")]
    eb = pd.read_csv(EBNER).set_index("protocol_id")
    pkt, meta, zones, outs = {}, dict(protos=PROTOS, traces={}, port_fail=[]), [], []
    for k, q in PROTOS.items():
        pairs = sorted(set(pp[pp.target == q].pair))
        want = ch["pair"] if ch["pair"] in pairs and T.has_rec(PATH, ch["pair"], q) else None
        rat = {}
        for p in pairs:
            if not T.has_rec(PATH, p, q):
                continue
            spv, *_ = T.pair_spk(T.FJ("CHR"), PATH, p, [q])
            r = T.load_rec(PATH, p, q)
            n, Tn = r["effcai"].shape
            t = np.asarray(r["t"], np.float64); dt = np.diff(t) / 1000.0
            args, info = T.prep(T.FJ("CHR"), PATH, r, spv)
            ref = T.ref_run(args, n)
            rat[p] = float(np.nanmean(T.ratio_hook(r, args, info, ref)["ratio_ctl"]))
            TR, EV, NEV, tr = T.trace_run(args, n, Tn, list(range(n)), info["cnt"])
            if not T.check(f"{k} {p}", ref, tr, TR, np.arange(n))["ok"]:
                meta["port_fail"].append(f"{k} {p}")
            C = TR[:, T.F_["c"]]
            for i in range(n):
                c = C[i, :-1]; td, tp = float(info["td"][i]), float(info["tp"][i])
                zones.append(dict(proto=k, pair=p, syn=int(r["syn"][i]), rho0=float(r["rho0"][i]),
                                  rho_end=float(ref["rho"][i]), td=td, tp=tp, cmax=float(C[i].max()),
                                  t_dep=float(dt[(c >= td) & (c < tp)].sum()), t_pot=float(dt[c >= tp].sum())))
            if want is None:
                want = p
            if p == want and k not in meta["traces"]:
                syns = np.asarray(r["syn"])
                hit = np.flatnonzero(syns == int(ch["syn"]))
                i = int(hit[0]) if len(hit) else int(np.argmin(np.abs(r["rho0"] - 0.0)))
                pre = np.sort(np.asarray(r["prespikes"], float).ravel())
                pkt[f"{k}__t"] = t - pre[0]; pkt[f"{k}__c"] = C[i].astype(np.float32)
                pkt[f"{k}__pre"] = pre - pre[0]
                pkt[f"{k}__post"] = np.sort(np.asarray(r["postspikes"], float).ravel()) - pre[0]
                meta["traces"][k] = dict(pair=p, syn=int(r["syn"][i]), td=float(info["td"][i]),
                                         tp=float(info["tp"][i]), rho0=float(r["rho0"][i]))
            del TR, EV
            T._REC.pop((PATH, p, q), None); del r; gc.collect(); T.drop_recs(); gc.collect()
        x = pd.Series(rat).dropna()
        outs.append(dict(proto=k, target=q, mean=x.mean(), sem=x.std(ddof=1) / np.sqrt(len(x)), n=len(x),
                         data_mean=eb.loc[q, "mean_ratio"], data_sem=eb.loc[q, "sem"],
                         source=f"{eb.loc[q, 'source_paper']}, {eb.loc[q, 'figure_or_table']}"))
        z = pd.DataFrame([w for w in zones if w["proto"] == k])
        print(f"{k} {q}: ratio {x.mean():.3f} (n {len(x)}; data {outs[-1]['data_mean']}) | per syn: t_dep "
              f"{z.t_dep.mean():.3f} s, t_pot {z.t_pot.mean():.3f} s, cmax/tp {(z.cmax / z.tp).median():.2f} "
              f"({time.time() - t0:.0f} s)", flush=True)
    pd.DataFrame(zones).to_csv(os.path.join(WORK, f"{TAG}_zones.csv"), index=False)
    pd.DataFrame(outs).to_csv(os.path.join(WORK, f"{TAG}_outcomes.csv"), index=False)
    np.savez_compressed(os.path.join(WORK, f"{TAG}.npz"), **pkt)
    json.dump(meta, open(os.path.join(WORK, f"{TAG}_meta.json"), "w"), indent=1, default=str)
    print(f"wrote figV5 ({time.time() - t0:.0f} s); port failures {meta['port_fail'] or 'none'}", flush=True)


if __name__ == "__main__":
    main()
