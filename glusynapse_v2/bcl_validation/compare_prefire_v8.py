"""Prefire BCL validation of a v8 fit: live GluSynapseV8 (live_v8.py) vs the offline v8 rule. compare_prefire_v7.py is
reused unchanged (imported, its offline_v7 / ecb_cols patched); the offline reference is v8_rec, the CPU port of the
gpu_v8_rho kernel (gate_src 2, gate_win > 0, v5_mode 2) on top of compare_prefire_v7.v7_rec, with its arithmetic:
  shaft  X = rint(1e3 x (1e3 x (shaft_cai - rest))) nM clipped to 0..65535, read as X x 1e-3 uM; rest = gpu_v8_rho
         _shaft_rest (copied: shaft cai at the last grid sample before min(pre, post spikes) - 1 ms);
  licence at sample k: lic = [t_k - t_last <= gate_win], t_last = last sample j < k with X_j > theta_G (fit v8 theta_G);
         pot = 0 when theta_G > 0 and not lic (then dep (1 - pot) gamma_d acts); the Vg licence is not used;
  eCB / veto as v7_rec; t_exact 1 (fit SET, or env V8_T_EXACT=0|1 to override) adds the gpu_v7x_rho exact-arrival
         corrections (W(t_a) = W - f ip, Q(t_a) = Q + f up, f = (t_k - t_a) / h_j of the earliest arrival).
gate_src 0 fits fall through to compare_prefire_v7 unchanged. Per record it adds: nsh_live / nsh_off (licence openings,
sums over the common synapses), nsh_syn_mismatch, tsh_live / tsh_off (licence open time, ms), tsh_maxdiff,
rest_maxdiff_nM (live cai_rest vs offline rest). --show-syn also prints the per-synapse gate table.

    python glusynapse_v2/bcl_validation/compare_prefire_v8.py --fit <json> --results <jsonl> --save <prefix> [--show-syn]
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_prefire_v7 as C7                           # noqa: E402
from compare_prefire_v5 import C, MV, rho_v4, batch_v2, PATHS   # noqa: E402
import live_v7                                            # noqa: E402
from numba import njit                                    # noqa: E402

_offline_v7, _ecb_cols = C7.offline_v7, C7.ecb_cols
SHOW = "--show-syn" in sys.argv


@njit(cache=False)
def v8_rec(E, S, X, F, hs, cnt, td, tp, rho0, gd, gp, rs, kt, thG, win, tex, isc, tauE1, dmin, Ae, dp0, thE, tauD, Q,
           veto, out_rho, out_d, out_necb, out_nveto, out_nsh, out_tsh):
    """compare_prefire_v7.v7_rec with the gpu_v8_rho G1 licence (X (n, T) shaft dCa uM as the kernel reads it) and the
    optional gpu_v7x_rho corrections (tex; F (n, T) arrival fraction f at the arrival samples)."""
    n, T = E.shape
    for i in range(n):
        a_td = td[i]; a_tp = tp[i]
        r = rho0[i]; d0 = dp0
        V = 0.0; chV = -1.0; aV = 1.0; bV = 0.0
        W = 0.0; Bg = 0.0; chD = -1.0; aD = 1.0
        ip = 0.0; up = 0.0
        clk = 0.0; tl = -1e30; lp = False; nsh = 0.0; tsh = 0.0
        ne = 0.0; nv = 0.0
        thEi = thE[i]
        for k in range(T):
            e = E[i, k] * 1.0
            s = S[i, k] * 1.0
            h = hs[k]
            hm = h * kt
            c = cnt[i, k]
            fx = F[i, k] if tex else 0.0
            pot = 1.0 if e > a_tp else 0.0
            dep = 1.0 if e > a_td else 0.0
            lic = (clk - tl) <= win
            if lic and not lp:
                nsh += 1.0
            if lic:
                tsh += hm
            lp = lic
            if pot > 0.0 and thG > 0.0 and not lic:
                pot = 0.0
            r = r + h * (-r * (1 - r) * (rs - r) + pot * gp * (1 - r) - dep * (1 - pot) * gd * r)
            if r < 0.0:
                r = 0.0
            elif r > 1.0:
                r = 1.0
            if c > 0.0:
                Bg = Bg + c
            trig = (W - fx * ip) > thEi
            if c > 0.0 and Ae > 0.0 and trig:
                if veto and (Q[i, k] + fx * up) > thEi:
                    nv += c
                else:
                    d0 = dmin + (d0 - dmin) * (1.0 - Ae) ** c
                    ne += c
            if hm != chV:
                chV = hm
                aV = np.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / isc
            wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
            ip = bV * (wb * s)
            W = aV * W + ip
            if hm != chD:
                chD = hm
                aD = np.exp(-hm / tauD)
            Bg = Bg * aD
            up = bV * s
            V = aV * V + up
            if X[i, k] > thG:
                tl = clk
            clk += hm
        out_rho[i] = r
        out_d[i] = d0
        out_necb[i] = ne
        out_nveto[i] = nv
        out_nsh[i] = nsh
        out_tsh[i] = tsh


def shaft_rest(r, sh):
    """gpu_v8_rho._shaft_rest (copied, the GPU module is not imported on CPU nodes)."""
    t = np.asarray(r["t"], np.float64)
    ev = np.concatenate([np.asarray(r["prespikes"], float).ravel(), np.asarray(r["postspikes"], float).ravel()])
    k = 0
    if ev.size:
        k = max(int(np.searchsorted(t, ev.min() - 1.0)) - 1, 0)
    return sh[:, k].astype(np.float64)


def arrival_F(r, n, T):
    """gpu_v7x_rho.arrival_frac for one record as an (n, T) array (f at the arrival samples, 0 elsewhere)."""
    t = np.asarray(r["t"], np.float64)
    arr = np.asarray(r["arr"]).reshape(n, -1)
    delay = MV.edge_params(r["syn"])["delay"] if "delay" not in r else r["delay"]
    ta = np.asarray(r["prespikes"], np.float64).ravel()[None, :] + np.asarray(delay, np.float64)[:, None]
    F = np.zeros((n, T))
    for ii in range(n):
        k = arr[ii]; ks = np.unique(k)
        first = np.full(T, np.inf); np.minimum.at(first, k, ta[ii])
        tk = t[ks]; hj = tk - t[np.maximum(ks - 1, 0)]
        f = np.where(hj > 0, (tk - first[ks]) / np.where(hj > 0, hj, 1.0), 0.0)
        F[ii, ks] = np.clip(f, 0.0, 1.0)
    return F


def offline_v8(path, fit, P, pairs, protos, conds):
    src = int(P.get("gate_src", 0))
    if src == 0:
        return _offline_v7(path, fit, P, pairs, protos, conds)
    assert src == 2 and float(P.get("gate_win", 100.0)) > 0, "compare v8: gate_src 2 with gate_win > 0 (G1) only"
    win = float(P.get("gate_win", 100.0)); thG = float(fit["v8"]["theta_G"])
    tex = bool(int(os.environ.get("V8_T_EXACT", P.get("t_exact", 0))))
    batch_v2.BASIS_DIR = C.BASIS[path]
    B = batch_v2.BatchV2(PATHS[path]["dirs"], protocols=sorted(protos), pairs=set(pairs), fast=False,
                         signals=("vdcc", "shaft_cai"))
    assert int(P["v5_mode"]) == 2, "v8 = v5_mode 2"
    gamma, _, sigma = rho_v4.opts(P); assert gamma == 1.0
    gd, gp = rho_v4.rates(P); thV = rho_v4.vamp(P)[1]
    kt = 1e3 * float(batch_v2.TAU_IND_GB); rs = float(batch_v2.RHO_STAR_GB)
    Tv = float(P.get("veto_T", 0.0)); tauE1 = float(P["tau_E1"]); isc = float(P["i_scale"])
    assert Tv > 0 or not tex, "t_exact needs veto_T > 0"
    thE0 = float(P["theta_eCB"]); um = live_v7.ue_map(fit, path)
    print(f"v8 offline {path}: gate_src 2, win {win} ms, theta_G {thG} uM, t_exact {int(tex)}", flush=True)
    out = {}
    for r in B.recs:
        n, T = r["effcai"].shape
        hs = np.zeros(T); hs[:T - 1] = np.diff(r["t"]) / 1000.0 / batch_v2.TAU_IND_GB
        cnt = np.zeros((n, T)); np.add.at(cnt, (np.broadcast_to(np.arange(n)[:, None], r["arr"].shape), r["arr"]), 1.0)
        td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], 1.0)
        pre, post = map(int, r["pair"].split("-"))
        uE = np.ones(n) if um is None else np.array([um[(pre, post, int(s))] for s in r["syn"]])
        thE = uE * thE0 if thE0 > 0 else np.full(n, thV)
        Q = C7.veto_Q(r, cnt, Tv, tauE1, isc) if Tv > 0 else np.zeros((n, T))
        assert r.get("shaft_cai") is not None, f"{r['pair']}/{r['proto']}: no shaft_cai"
        sh = np.asarray(r["shaft_cai"], np.float64)
        rest = shaft_rest(r, sh)
        du = (sh - rest[:, None]) * 1e3                                    # as GPUModelV8.__init__
        X = np.clip(np.rint(du * 1e3), 0, 65535).astype(np.uint16).astype(np.float64) * 1e-3
        F = arrival_F(r, n, T) if tex else np.zeros((1, 1))
        rho_f = np.empty(n); d_f = np.empty(n); ne = np.empty(n); nv = np.empty(n); nsh = np.empty(n); tsh = np.empty(n)
        v8_rec(np.ascontiguousarray(r["effcai"], dtype=np.float32), np.ascontiguousarray(r["vdcc"], dtype=np.float32),
               X, F, hs, cnt, np.asarray(td, np.float64), np.asarray(tp, np.float64), r["rho0"].astype(np.float64), gd,
               gp, rs, kt, thG, win, tex, isc, tauE1, float(P["dpre_min"]), float(P["A_eCB"]), 0.0 + float(P["dpre0"]),
               np.ascontiguousarray(thE, np.float64), float(P["tau_d_NMDA"]), Q, Tv > 0, rho_f, d_f, ne, nv, nsh, tsh)
        b = B.basis(r)
        ncev = np.isfinite(r["cev"]).sum(axis=1) if r.get("cev") is not None else None
        for c in conds:
            frozen = c in ("post_nmdar", "nmdar_block")
            blocked = c in ("mglu_block", "nmdar_block")
            rc = r["rho0"] if frozen else rho_f
            d = (np.zeros(n) if c == "nmdar_block" else np.full(n, 0.0 + float(P["dpre0"])) if c == "mglu_block"
                 else d_f)
            out[(r["pair"], r["proto"], c)] = dict(
                syn=np.asarray(r["syn"]), rho0=r["rho0"], rho_obs=np.asarray(r["rho_obs"], float), rho=np.asarray(rc),
                dpre=np.asarray(d), ncev=ncev, b=b, ratio=rho_v4.ratio(b, r["rho0"], rc, d, 0.0 if frozen else sigma),
                necb=np.zeros(n) if blocked else ne, nveto=np.zeros(n) if blocked else nv, uE=uE,
                nsh=nsh, tsh=tsh, rest=rest)
        r["shaft_cai"] = None
    return out


def ecb_cols(x, o):
    """compare_prefire_v7.ecb_cols + the gate columns (and the per-synapse gate table with --show-syn)."""
    d = _ecb_cols(x, o)
    if not d or "v8" not in x or "nsh" not in o:
        return d
    live = {int(s): i for i, s in enumerate(x["syn"])}
    idx = np.array([live.get(int(s), -1) for s in o["syn"]]); ok = idx >= 0
    g = x["v8"]
    nl = np.array(g["nsh"])[idx[ok]]; tl = np.array(g["tsh"])[idx[ok]]; rl = np.array(g["cai_rest"])[idx[ok]]
    no = o["nsh"][ok]; to = o["tsh"][ok]; ro = o["rest"][ok]
    d.update(nsh_live=float(nl.sum()), nsh_off=float(no.sum()), nsh_syn_mismatch=int(np.sum(nl != no)),
             tsh_live=float(tl.sum()), tsh_off=float(to.sum()), tsh_maxdiff=float(np.max(np.abs(tl - to))),
             rest_maxdiff_nM=float(1e6 * np.max(np.abs(rl - ro))))
    if SHOW:
        print(f"\n{x['path']} {x['pair']} {x['proto']} {x['cond']} gate: syn | rest nM live/off | nsh live/off | "
              f"tsh ms live/off | rho live/off (theta_G {g['theta_sh_uM']:.4g} uM, t_rest {g['t_rest']:.2f} ms)")
        st = x["end"]
        for j, s in enumerate(o["syn"]):
            i = live.get(int(s))
            if i is None:
                continue
            print(f"  {int(s)} | {1e6 * g['cai_rest'][i]:.2f}/{1e6 * o['rest'][j]:.2f} | {g['nsh'][i]:.0f}/{o['nsh'][j]:.0f}"
                  f" | {g['tsh'][i]:.1f}/{o['tsh'][j]:.1f} | {st['rho_GB'][i]:.4f}/{o['rho'][j]:.4f}")
    return d


C7.offline_v7 = offline_v8
C7.ecb_cols = ecb_cols


def main():
    C7.main()
    a = sys.argv
    f = a[a.index("--save") + 1] + "_records.csv"
    if os.path.isfile(f):
        import pandas as pd
        df = pd.read_csv(f)
        if "nsh_live" in df:
            d = df[df.ok.fillna(False).astype(bool)]
            print(f"\nv8 gate: licence openings live {d.nsh_live.sum():.0f} vs offline {d.nsh_off.sum():.0f} "
                  f"({int(d.nsh_syn_mismatch.sum())} syn differ), open time live {d.tsh_live.sum():.1f} vs offline "
                  f"{d.tsh_off.sum():.1f} ms (max per-syn diff {d.tsh_maxdiff.max():.2f} ms), rest max diff "
                  f"{d.rest_maxdiff_nM.max():.3f} nM")
            for p, g in d.groupby("path"):
                print(f"v8 gate {p}: {len(g)} records, openings live {g.nsh_live.sum():.0f} vs offline "
                      f"{g.nsh_off.sum():.0f} ({int(g.nsh_syn_mismatch.sum())} syn differ), open time live "
                      f"{g.tsh_live.sum():.1f} vs offline {g.tsh_off.sum():.1f} ms (max per-syn diff "
                      f"{g.tsh_maxdiff.max():.2f} ms), rest max diff {g.rest_maxdiff_nM.max():.3f} nM, rho disagreements "
                      f"{int(g.rho_disagree.sum())} / {int(g.n_common.sum())} syn")
            if len(d) <= 40:
                with pd.option_context("display.width", 250, "display.max_columns", 30):
                    print(d[["path", "pair", "proto", "cond", "ratio_basis_live", "ratio_off", "rho_disagree", "nsh_live",
                             "nsh_off", "tsh_live", "tsh_off", "tsh_maxdiff", "rest_maxdiff_nM"]].round(4)
                          .to_string(index=False))
    f = a[a.index("--save") + 1] + "_targets.csv"
    if os.path.isfile(f):
        import pandas as pd
        o = pd.read_csv(f)
        for p, g in o.groupby("path"):
            fi = g[~g.validation_only.astype(bool)]
            print(f"v8 chi2 {p}: {len(fi)} fitted targets, BCL live {fi.z_live2.sum():.2f}, offline same records "
                  f"{fi.z_off2.sum():.2f}, fit csv {np.nansum(fi.z_fit2):.2f}; PASS {int(g.PASS.sum())}/{len(g)} "
                  f"(incl. {int(g.validation_only.astype(bool).sum())} validation-only)")


if __name__ == "__main__":
    main()
