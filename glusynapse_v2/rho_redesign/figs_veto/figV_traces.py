"""Traces for figV1 (Chindemi refit) and figV2 (B4, veto_peak_k 0.75): why the eCB step needs the timed veto.

Uses tv12_lat.py = a verbatim copy of /scratch/dhuruva/latest_figs/traces/traces_v12_lat.py (the figLC / figLD trace
tool: v12_tr = compare_prefire_v12.v12_rec + recording; port check against v12_rec on every traced record).
Nothing in the rule or the parameters is changed. Three parameter sets:
  B4  = /scratch/dhuruva/perpair/w_z14b4k75_u97_fit.json (4V + B4, veto_peak_k 0.75, d_NO_max 0.14; kernel gpu_v12_rho)
  KOV = B4 with veto_peak_k 1e30 (veto off, no refit; as latest_figs/variants/KO_veto: C never fires either)
  CHR = /scratch/dhuruva/latest_figs/variants/chr/CHR_s5.json (Chindemi form refit, the orange model of figL0)
Protocols (L5->L5, s2g0321 records): sjostrom_0.1hz_dt-10ms, 10Hz_-10ms (Markram), 10Hz_10ms (Markram).

1. Scan: every synapse of every pair of the u97 per-pair table (control rows of the three targets) through v12_rec
   under B4 / KOV / CHR; per synapse end state, per pair EPSP ratio (tv12_lat.ratio_hook = rho_v4.ratio on the pair
   basis, control lane). -> WORK/figV_pop.csv
2. Outcomes: per (protocol, model) mean +/- SEM over pairs of the CPU ratio, next to the GPU val csv (pred, pred_sem)
   and the GPU per-pair mean. -> WORK/figV_outcomes.csv
3. Synapse choice (one synapse for all columns and both figures): synapses present in all three records that step at
   0.1 Hz -10 and at 10 Hz -10 under B4 and have a vetoed trigger at 10 Hz +10 under B4; contrast
   |d_eCB_end(KOV) - d_eCB_end(B4)| at 10 Hz +10; the one closest to the median contrast among contrast > 0.01
   (the figLD rule). Fallbacks drop the -10 conditions one by one (logged).
4. Traces of that synapse, 3 protocols x 3 models -> WORK/figV_traces.npz (+ figV_meta.json). Spine Ca (cai) is not in
   the records; it is recovered from effcai: the extractor computed effcai' = -effcai/tau + (cai - min_ca) exactly, so
   the interval mean cai = (E[k+1] - A E[k]) / (tau (1 - A)) + min_ca, A = exp(-dt_k / tau).

  python glusynapse_v2/rho_redesign/figs_veto/figV_traces.py
MEASURED 22479950: 2:21, 3.60 GB of 3700M (97%), CPU eff 50% -> next run 4500M, 0:15.
(env as traces_lat_LCLD: BCL_EMODEL, ANALYTICAL_BASIS_DIR, L23_BASIS_DIR, L23L23_BASIS_DIR, CEXP_DIR)
"""
import copy, gc, json, os, sys, time, traceback
import argparse
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tv12_lat as T                                      # noqa: E402

import importlib.util                                     # noqa: E402
_sp = importlib.util.spec_from_file_location(
    "am_constants", "/lustre09/project/6070394/dhuruva/plastyfire/analytical_method/constants.py")
_am = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(_am)
MIN_CA_CR, TAU_EFFCA_GB = float(_am.MIN_CA_CR), float(_am.TAU_EFFCA_GB)   # 70e-6 mM, 278.3177658387 ms

PP = "/scratch/dhuruva/perpair"
B4J = f"{PP}/w_z14b4k75_u97_fit.json"
CHRJ = "/scratch/dhuruva/latest_figs/variants/chr/CHR_s5.json"
GPUVAL = dict(B4=(f"{PP}/w_z14b4k75_u97_val.csv", f"{PP}/w_z14b4k75_u97_val_perpair.csv"),
              CHR=("/scratch/dhuruva/latest_figs/variants/CHR_s5_val.csv",
                   "/scratch/dhuruva/latest_figs/variants/CHR_s5_val_perpair.csv"))
PROTOS = dict(s01="sjostrom_0.1hz_dt-10ms", m10="10Hz_-10ms", p10="10Hz_10ms")
FITS = ("B4", "KOV", "CHR")
PATH = "L5"
TRK = ("c", "vdcc", "pool", "d_eCB", "d_NO", "rho", "clk")


def kov(fit):
    f = copy.deepcopy(fit)
    s = json.loads(f["args"]["set"]); s["veto_peak_k"] = 1e30; f["args"]["set"] = json.dumps(s)
    f["pre"]["veto_peak_k"] = 1e30
    for b in ("v10", "v11"):
        if isinstance(f.get(b), dict) and isinstance(f[b].get("params"), dict) and "veto_peak_k" in f[b]["params"]:
            f[b]["params"]["veto_peak_k"] = 1e30
    f["variant"] = dict(name="KOV", base=B4J, set_changes={"veto_peak_k": 1e30}, refit=False)
    return f


def cai_from_effcai(E, t):
    """Interval-mean spine Ca (uM) from effcai on the record grid (exact inverse of extract.effcai_from_cai)."""
    E = np.asarray(E, np.float64); dt = np.diff(np.asarray(t, np.float64))
    A = np.exp(-dt / TAU_EFFCA_GB)
    u = (E[1:] - A * E[:-1]) / (TAU_EFFCA_GB * (1.0 - A))
    return np.append((u + MIN_CA_CR) * 1e3, np.nan)          # mM -> uM; value of step k = mean over [t_k, t_k+1)


def choose(w, key, mask, eps, what):
    """tv12_lat.choose on a (pair, syn)-keyed table: among mask & contrast > eps, the one closest to the median."""
    c = w[mask & (w[key] > eps)]
    info = dict(contrast=what, eps=eps, n_syn=int(len(w)), n_cand=int(len(c)),
                rule="closest to the median contrast among synapses where the veto acts (contrast > eps)")
    if c.empty:
        return None, info
    med = float(c[key].median())
    row = c.iloc[int(np.argmin(np.abs(c[key].to_numpy() - med)))]
    info.update(median=med, chosen=float(row[key]), max=float(c[key].max()), q25=float(c[key].quantile(.25)),
                q75=float(c[key].quantile(.75)))
    return row, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/scratch/dhuruva/figs_veto")
    A = ap.parse_args()
    os.makedirs(A.out, exist_ok=True)
    T.A = argparse.Namespace(out=A.out); T.SEL, T.STATS, T.CHOICE, T.VFILES = [], {}, {}, {}
    b4 = json.load(open(B4J))
    T._FJ.update(B4=b4, KOV=kov(b4), CHR=json.load(open(CHRJ)))
    t0 = time.time()

    # ---- pairs: control rows of the three targets in the u97 per-pair table
    pp = pd.read_csv(GPUVAL["B4"][1])
    pp = pp[(pp.pathway == "l5") & (pp.condition == "control")]
    pairs = {k: set(pp[pp.target == q].pair) for k, q in PROTOS.items()}
    allp = sorted(set().union(*pairs.values()))
    print(f"pairs per protocol { {k: len(v) for k, v in pairs.items()} }, union {len(allp)}", flush=True)

    # ---- 1. scan
    rows, errs = [], set()
    for ip, p in enumerate(allp):
        have = [k for k in PROTOS if p in pairs[k] and T.has_rec(PATH, p, PROTOS[k])]
        if not have:
            continue
        spv, tiers, *_ = T.pair_spk(T.FJ("B4"), PATH, p, [PROTOS[k] for k in have])
        for k in have:
            q = PROTOS[k]
            r = T.load_rec(PATH, p, q)
            n = r["effcai"].shape[0]
            for fk in FITS:
                try:
                    args, info = T.prep(T.FJ(fk), PATH, r, spv)
                    o = T.ref_run(args, n)
                    dc, _ = T.lanes(info, o["d"], o["dn"])
                    rr = T.ratio_hook(r, args, info, o)
                except Exception as e:
                    if fk not in errs:
                        traceback.print_exc()
                    errs.add(fk)
                    print(f"  {p} {q} {fk}: ERROR {type(e).__name__}: {e}", flush=True)
                    continue
                for i in range(n):
                    rows.append(dict(pair=p, proto=k, fit=fk, i=i, syn=int(r["syn"][i]), rho0=float(r["rho0"][i]),
                                     dctl=float(dc[i]), ratio=float(rr["ratio_ctl"][i]),
                                     **{kk: float(o[kk][i]) for kk in T.OUTK}))
            T._REC.pop((PATH, p, q), None); del r; gc.collect()
        T.drop_recs(); gc.collect()
        print(f"  scan {ip + 1}/{len(allp)} {p}: {len(rows)} rows ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(A.out, "figV_pop.csv"), index=False)

    # ---- 2. outcomes
    out = []
    for k, q in PROTOS.items():
        for fk in FITS:
            x = df[(df.proto == k) & (df.fit == fk)].groupby("pair").ratio.first().dropna()
            row = dict(proto=k, target=q, fit=fk, cpu_mean=x.mean(), cpu_sem=x.std(ddof=1) / np.sqrt(len(x)),
                       cpu_n=len(x))
            if fk in GPUVAL:
                v = pd.read_csv(GPUVAL[fk][0]); v = v[(v.target == q) & (v.condition == "control")].iloc[0]
                g = pd.read_csv(GPUVAL[fk][1])
                g = g[(g.pathway == "l5") & (g.target == q) & (g.condition == "control")].set_index("pair").ratio
                cm = pd.concat([x.rename("cpu"), g.rename("gpu")], axis=1, join="inner")
                row.update(gpu_pred=v.pred, gpu_sem=v.pred_sem, gpu_n=int(v.n_pairs), data_mean=v.target_mean,
                           data_sem=v.target_sem, gpu_pp_mean=g.mean(), gpu_pp_sem=g.std(ddof=1) / np.sqrt(len(g)),
                           cpu_gpu_maxdiff=float((cm.cpu - cm.gpu).abs().max()) if len(cm) else np.nan,
                           cpu_gpu_n=len(cm))
            out.append(row)
    oc = pd.DataFrame(out)
    oc.to_csv(os.path.join(A.out, "figV_outcomes.csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print("\nOUTCOMES (CPU port pair ratios vs GPU val)\n" + oc.to_string(index=False, float_format="%.4f"),
              flush=True)

    # ---- 3. synapse choice
    def sel(k, fk, tag):
        x = df[(df.proto == k) & (df.fit == fk)].set_index(["pair", "syn"])
        return x.drop(columns=["proto", "fit"]).add_suffix("_" + tag)
    w = sel("s01", "B4", "s").join(sel("m10", "B4", "m"), how="inner").join(sel("p10", "B4", "on"), how="inner") \
        .join(sel("p10", "KOV", "off"), how="inner").reset_index()
    w["dd"] = (w.d_off - w.d_on).abs()
    print(f"\nsynapses in all three records: {len(w)} ({w.pair.nunique()} pairs); vetoed at 10Hz +10 under B4: "
          f"{int((w.nv_on > 0).sum())}; step at 0.1Hz -10: {int((w.ne_s > 0).sum())}; step at 10Hz -10: "
          f"{int((w.ne_m > 0).sum())}; |d_KOV - d_B4| > 0.01 at 10Hz +10: {int((w.dd > 0.01).sum())}", flush=True)
    masks = [("steps at 0.1Hz -10 and at 10Hz -10, vetoed trigger at 10Hz +10", (w.ne_s > 0) & (w.ne_m > 0) & (w.nv_on > 0)),
             ("steps at 10Hz -10, vetoed trigger at 10Hz +10", (w.ne_m > 0) & (w.nv_on > 0)),
             ("vetoed trigger at 10Hz +10", w.nv_on > 0)]
    row = None
    for desc, m in masks:
        row, ci = choose(w, "dd", m, 0.01, "|d_eCB_end(veto off) - d_eCB_end(B4)| at 10Hz +10")
        if row is not None:
            ci["mask"] = desc
            break
        print(f"  no candidate with mask '{desc}'", flush=True)
    assert row is not None, "no synapse found"
    pr, syn = row["pair"], int(row["syn"])
    ci.update(pair=pr, syn=syn, rho0=float(row["rho0_s"]))
    print(f"CHOSEN {pr} syn {syn}: {json.dumps(ci, default=str)}", flush=True)
    print("  chosen synapse end states:\n" + w[(w.pair == pr) & (w.syn == syn)].T.to_string(), flush=True)

    # ---- 4. traces
    spv, tiers, *_ = T.pair_spk(T.FJ("B4"), PATH, pr, list(PROTOS.values()))
    pkt, meta = {}, dict(choice=ci, fits=dict(B4=B4J, KOV="B4 with veto_peak_k 1e30", CHR=CHRJ), protos=PROTOS,
                         traces={}, port_checks={})
    for k, q in PROTOS.items():
        r = T.load_rec(PATH, pr, q)
        i = int(np.flatnonzero(np.asarray(r["syn"]) == syn)[0])
        n, Tn = r["effcai"].shape
        t = np.asarray(r["t"], np.float64)
        delay = float(np.asarray(T.MV.edge_params(r["syn"])["delay"], float)[i])
        pre = np.sort(np.asarray(r["prespikes"], float).ravel()); post = np.sort(np.asarray(r["postspikes"], float).ravel())
        pkt[f"{k}__t"] = t; pkt[f"{k}__pre"] = pre; pkt[f"{k}__post"] = post
        pkt[f"{k}__cai"] = cai_from_effcai(r["effcai"][i], t).astype(np.float32)
        cr = pkt[f"{k}__cai"]
        print(f"\n{k} {q}: {n} synapses, {Tn} steps, t {t[0]:.1f}..{t[-1]:.1f} ms, {len(pre)} pre / {len(post)} post "
              f"spikes, delay {delay:.2f} ms; pre[:6] {np.round(pre[:6], 1).tolist()} post[:8] "
              f"{np.round(post[:8], 1).tolist()}; recovered cai q01/50/max {np.nanquantile(cr, .01):.4f} / "
              f"{np.nanmedian(cr):.4f} / {np.nanmax(cr):.3f} uM", flush=True)
        for fk in FITS:
            args, info = T.prep(T.FJ(fk), PATH, r, spv)
            ref = T.ref_run(args, n)
            TR, EV, NEV, tr = T.trace_run(args, n, Tn, [i], info["cnt"])
            chk = T.check(f"{k} {fk}", ref, tr, TR, np.asarray([i]))
            meta["port_checks"][f"{k}__{fk}"] = chk
            dctl, _ = T.lanes(info, TR[0, T.F_["d_eCB"]], TR[0, T.F_["d_NO"]])
            for f in TRK:
                pkt[f"{k}__{fk}__{f}"] = TR[0, T.F_[f]].astype(np.float64 if f in ("clk", "rho", "d_eCB", "d_NO")
                                                                   else np.float32)
            pkt[f"{k}__{fk}__dctl"] = np.asarray(dctl, np.float64)
            ne = int(NEV[0]); pkt[f"{k}__{fk}__EV"] = EV[0, :ne]
            # peak own VDCC current in the veto window (t_a + vt0, t_a + vtv] per arrival (diagnostic; kernel flag rules)
            clk = TR[0, T.F_["clk"]]; S = TR[0, T.F_["vdcc"]]
            h = np.r_[np.diff(clk), 0.0]
            ta = EV[0, :ne, T.EV_FIELDS.index("t_a")]
            qa = np.array([np.nanmax(np.where((clk + h > a + info["vt0"]) & (clk < a + info["vtv"]) & (h > 0), S, -np.inf))
                           for a in ta]) if ne else np.zeros(0)
            pkt[f"{k}__{fk}__qa"] = qa
            sc = dict(td=float(info["td"][i]), tp=float(info["tp"][i]), thE=float(info["thE"][i]),
                      thP=float(info["thP"][i]), uE=float(info["uE"][i]), spk=float(info["spk"][i]),
                      dmin=info["dmin"], dnmax=float(info["dnmax"]) if info["dnmax"] < 1e29 else None,
                      vt0=info["vt0"], vtv=info["vtv"], k_peak=info["k_peak"], theta_eCB=info["theta_eCB"],
                      nm=int(info["nm"]), Ae=info["Ae"], rho0=float(r["rho0"][i]), delay=delay, i=i,
                      end={kk: float(ref[kk][i]) for kk in T.OUTK}, dctl_end=float(dctl[-1]))
            meta["traces"][f"{k}__{fk}"] = sc
            E = pkt[f"{k}__{fk}__EV"]; ix = {f: T.EV_FIELDS.index(f) for f in T.EV_FIELDS}
            print(f"  {fk}: rho {sc['rho0']:.0f} -> {sc['end']['rho']:.3f}, d_eCB end {sc['end']['d']:.3f}, d_NO end "
                  f"{sc['end']['dn']:.3f}, control d end {sc['dctl_end']:.3f}; arrivals {ne}, triggers "
                  f"{int(E[:, ix['trig']].sum())}, steps {int(E[:, ix['step']].sum())}, vetoed {int(E[:, ix['veto']].sum())}"
                  f", annulled {int(E[:, ix['annul']].sum())}; theta_d {sc['td']:.4g} theta_p {sc['tp']:.4g}; "
                  f"thE {sc['thE']:.4g} (uE {sc['uE']:.4g}), thP {sc['thP']:.4g}", flush=True)
            if fk == "B4":
                print("    j   t_a(ms rel pre0)  pool/P1  trig  veto  step  annul  C  qa/thP  d_eCB_after", flush=True)
                for j in range(min(ne, 12 if k != "s01" else 4)):
                    print(f"    {j:2d}  {E[j, ix['t_a']] + t[0] - pre[0]:9.1f}  {E[j, ix['trigval']] / sc['uE']:8.3f}  "
                          f"{int(E[j, ix['trig']]):4d}  {int(E[j, ix['veto']]):4d}  {int(E[j, ix['step']]):4d}  "
                          f"{int(E[j, ix['annul']]):5d}  {int(E[j, ix['C']]):1d}  {qa[j] / sc['thP']:6.3f}  "
                          f"{E[j, ix['d_eCB']]:8.3f}", flush=True)
            del TR
        T._REC.pop((PATH, pr, q), None); del r; gc.collect()
    meta["ev_fields"] = list(T.EV_FIELDS); meta["min_ca_mM"] = MIN_CA_CR; meta["tau_effca"] = TAU_EFFCA_GB
    np.savez_compressed(os.path.join(A.out, "figV_traces.npz"), **pkt)
    json.dump(meta, open(os.path.join(A.out, "figV_meta.json"), "w"), indent=1, default=str)
    print(f"\nwrote {A.out}/figV_traces.npz, figV_meta.json, figV_outcomes.csv, figV_pop.csv "
          f"({time.time() - t0:.0f} s); port failures {T.FAIL or 'none'}; scan errors {sorted(errs) or 'none'}", flush=True)
    if T.FAIL:
        sys.exit(3)


if __name__ == "__main__":
    main()
