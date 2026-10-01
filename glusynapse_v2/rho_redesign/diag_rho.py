"""Diagnosis of the post rho rule (RHO_REDESIGN.md step 1). CPU/numba, no jax. Fixed fit, no refit.

Per synapse (every record of every protocol in --dirs): c_pre, c_post, theta_d, theta_p, effcai peak, time above
theta_p and in the depression band (theta_d < effcai <= theta_p), rho0, rho_f (full ODE), dpre (control / eCB off /
NO off), soma path distance and kind where local_t/out has it (L5 subset24).
Per record: EPSP ratio for the decompositions
    full       rho_f, dpre control                    (= the fit's control prediction)
    rho_only   rho_f, dpre 0                          (post rule alone)
    pre_only   rho0,  dpre control                    (pre rule alone)
    no_ecb     rho_f, dpre with A_mglu 0              (= mglu_block)
    no_no      rho_f, dpre with A_NO 0                (= no_block)
    cont       rho_f continuous in the readout (no >= 0.5 binarisation), dpre control
Example traces (effcai / theta_p around the first pairings) for a few protocols go to <save>_traces.npz.

    ANALYTICAL_BASIS_DIR=... python glusynapse_v2/rho_redesign/diag_rho.py --fit <json> --dirs <d1,d2> --save <prefix> [--l23]
"""
import argparse, glob, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2)
import batch_v2                                  # noqa: E402
from batch_v2 import BatchV2, RHO_BINARY_THRESHOLD   # noqa: E402
import model_v2 as MV                            # noqa: E402

TRACE_PROTOS = ("sjostrom_0.1hz_dt-10ms", "sjostrom_0.1hz_dt+10ms", "sjostrom_20hz_dt-10ms", "sjostrom_20hz_dt+10ms",
                "sjostrom_50hz_dt-10ms", "sjostrom_50hz_dt+10ms", "letzkus_3ap_200hz_dt-10ms", "letzkus_3ap_200hz_dt+10ms",
                "letzkus_1ap_dt+10ms")


def distances():
    """syn id -> (path distance um, kind) from local_t/out (og variant, L5 subset24 synapses)."""
    out = {}
    for f in glob.glob(os.path.join(V2, "local_t", "out", "*.npz")):
        z = np.load(f)
        for s, d, k in zip(z["og|ap1|sid"], z["og|ap1|dist"], z["og|ap1|kind"]):
            out[int(s)] = (float(d), str(k))
    return out


def epsp_cont(b, rho, dpre):
    """PairBasis.epsp with rho continuous (linear in rho, as the basis is linear in the flipped state)."""
    rb = np.clip(np.asarray(rho, float), 0, 1)
    u0 = b.Use_d + rb * (b.Use_p - b.Use_d)
    scale = np.minimum(1.0, u0 * (1.0 + dpre)) / u0
    mean = float(np.sum((b.e0 * b.w + rb * b.delta) * scale))
    var = float(np.sum(rb * b.svar) + (1 - rb.sum()) ** 2 * b.v0)
    return mean, np.sqrt(max(var, 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--dirs", required=True)
    ap.add_argument("--save", required=True); ap.add_argument("--l23", action="store_true")
    ap.add_argument("--groups", default=None, help="restrict to the protocols of these target groups")
    a_ = ap.parse_args()
    fit = json.load(open(a_.fit)); fa = fit["args"]
    fil = json.loads(fa["filters"]); P = {**MV.DEFAULTS, **fil, **fit["pre"]}; A = fit["a"]
    pairs = None
    loc = {}
    if a_.l23:
        g = pd.read_csv(os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv"))
        g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
        pairs = set(g.pair)
        loc = {p: dict(letzkus_distal=bool(ld), sh_distal=bool(sd), path_mean=float(pm))
               for p, ld, sd, pm in zip(g.pair, g.letzkus_distal, g.sh_distal, g.path_mean)}
    protos = None
    if a_.groups:
        from targets import load_targets
        protos = sorted({k[0].split("@")[0] for k in load_targets(tuple(a_.groups.split(",")))})
    B = BatchV2(a_.dirs.split(","), protocols=protos, pairs=pairs, fast=False, signals=("vdcc",))
    th = [B.thetas(r, A) for r in B.recs]
    rho_all = B.rho_full_all(np.concatenate([t[0] for t in th]), np.concatenate([t[1] for t in th]))
    feats = B.features(P)
    st = np.load(os.path.join(V2, "syn_section_type.npz")); styp = dict(zip(st["syn"].tolist(), st["section_type"].tolist()))
    dist = distances()
    syn_rows, rec_rows, traces = [], [], {}
    for i, r in enumerate(B.recs):
        td, tp = th[i]; rho = rho_all[r["sl"]]; E = r["effcai"]; h = np.diff(r["t"]) / 1000.0
        above_p = ((E[:, :-1] > tp[:, None]) * h).sum(1)
        band = (((E[:, :-1] > td[:, None]) & (E[:, :-1] <= tp[:, None])) * h).sum(1)
        tT, K = feats[i]
        d = {k: MV.dpre_final(tT, K, am, ano, P["dpre_min"], P["dpre_max"], P["dpre0"])
             for k, am, ano in (("ctrl", P["A_mglu"], P["A_NO"]), ("noecb", 0.0, P["A_NO"]), ("nono", P["A_mglu"], 0.0))}
        b = B.basis(r); z = np.zeros(len(rho))
        rt = dict(full=b.ratio(r["rho0"], rho, d["ctrl"]), rho_only=b.ratio(r["rho0"], rho, z),
                  pre_only=b.ratio(r["rho0"], r["rho0"], d["ctrl"]), no_ecb=b.ratio(r["rho0"], rho, d["noecb"]),
                  no_no=b.ratio(r["rho0"], rho, d["nono"]))
        bm, bs = b.epsp(r["rho0"], z); am_, _ = epsp_cont(b, rho, d["ctrl"])
        rt["cont"] = (am_ / bm) * (1 + min((bs / bm) ** 2, 0.25)) if bm else np.nan
        L = loc.get(r["pair"], {})
        rec_rows.append(dict(pair=r["pair"], proto=r["proto"], n_syn=len(rho), **rt, **L))
        for k in range(len(rho)):
            s = int(r["syn"][k]); dd = dist.get(s, (np.nan, ""))
            syn_rows.append(dict(pair=r["pair"], proto=r["proto"], syn=s, c_pre=r["c_pre"][k], c_post=r["c_post"][k],
                                 td=td[k], tp=tp[k], peak=float(r["peak"][k]), t_above_tp=above_p[k], t_band=band[k],
                                 rho0=r["rho0"][k], rho_f=rho[k], dpre=d["ctrl"][k], dpre_noecb=d["noecb"][k],
                                 dpre_nono=d["nono"][k], dist=dd[0], kind=dd[1], sec_type=styp.get(s, -1), **L))
        if r["proto"] in TRACE_PROTOS and r["proto"] not in traces:
            t = r["t"]; t0 = float(r["prespikes"][0]) if len(r["prespikes"]) else float(t[0])
            m = (t >= t0 - 200) & (t <= t0 + 2500)
            traces[r["proto"]] = dict(pair=r["pair"], t=t[m] - t0, E=E[:, m], td=td, tp=tp, c_pre=r["c_pre"],
                                      c_post=r["c_post"], pre=np.asarray(r["prespikes"]) - t0,
                                      post=np.asarray(r["postspikes"]) - t0, rho_f=rho, rho0=r["rho0"])
    S = pd.DataFrame(syn_rows); R = pd.DataFrame(rec_rows)
    S.to_csv(a_.save + "_syn.csv", index=False); R.to_csv(a_.save + "_rec.csv", index=False)
    np.savez_compressed(a_.save + "_traces.npz", **{f"{p}|{k}": v for p, dct in traces.items() for k, v in dct.items()
                                                  if k != "pair"}, pairs=np.array([f"{p}|{d['pair']}" for p, d in traces.items()]))
    S["floor"] = S.c_post < 2e-4; S["m_p"] = S.peak / S.tp; S["m_d"] = S.peak / S.td
    S["up"] = (S.rho0 < 0.5) & (S.rho_f >= 0.5); S["down"] = (S.rho0 >= 0.5) & (S.rho_f < 0.5)
    pd.set_option("display.width", 250)
    print("=== per protocol: decomposition (mean over records)")
    print(R.groupby("proto")[["full", "rho_only", "pre_only", "no_ecb", "no_no", "cont"]].mean().round(3).to_string())
    print("=== per protocol x floor: margins and transitions")
    g = S.groupby(["proto", "floor"])
    T = pd.DataFrame(dict(n=g.size(), m_p_med=g.m_p.median(), m_p_q10=g.m_p.quantile(0.1), m_d_med=g.m_d.median(),
                          f_above_tp=g.m_p.apply(lambda x: (x > 1).mean()), t_tp_med=g.t_above_tp.median(),
                          t_band_med=g.t_band.median(),
                          P_up=g.apply(lambda x: x.up.sum() / max((x.rho0 < 0.5).sum(), 1)),
                          P_down=g.apply(lambda x: x.down.sum() / max((x.rho0 >= 0.5).sum(), 1)),
                          dpre=g.dpre.mean()))
    print(T.round(3).to_string())
    print(f"saved {a_.save}_syn.csv / _rec.csv / _traces.npz  ({len(S)} synapse rows, {len(R)} records)")


if __name__ == "__main__":
    main()
