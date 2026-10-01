"""dt 0 live / offline eCB mismatch of the v7 BCL validation (BCL_W2_N3.md). Per own arrival of the L5 records of
--protos it reads the offline trigger pool W and the veto charge Q at the kernel's arrival sample k (first grid sample
at or after t_a = prespike + delay) and re-scores three conventions against the live counts (merged jsonl):
  V0 kernel   trig W(t_k) > thE, veto Q[t_k, t_k + T] > thE (must reproduce compare_prefire_v7's necb_off / nveto_off);
  V1 exact    the live convention: W(t_a) = W(t_k) - (t_k - t_a)/h_j x (bin j = [t_k-1, t_k] W input), Q from t_a
              (+ the same fraction of bin j, unweighted);
  V2 prev     arrival on the sample before t_a: W(t_k-1), Q + all of bin j.
Prints per protocol: gap t_k - t_a, nearest post spike - t_a, W/thE at arrivals, triggered / step / veto totals per
convention vs live, and the number of synapses whose (necb, nveto) equal live.

    python glusynapse_v2/bcl_validation/diag_dt0_v7.py --fit <json> --results <jsonl> --protos a,b,...
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import compare_prefire_v7 as C7                              # noqa: E402
from compare_prefire_v7 import C, MV, rho_v4, batch_v2, PATHS, V5_DEFAULTS, live_v7   # noqa: E402
from numba import njit                                       # noqa: E402


@njit(cache=False)
def pools(S, hs, cnt, kt, isc, tauE1, tauD, Wk, Wp, inc, unw, hp):
    """W of the v7 kernel (weighted, v5_mode 2) and, at each arrival sample k: W before step k (Wk), W before step k-1
    (Wp), the weighted (inc) and unweighted (unw) input of step k-1, and its length hp (ms)."""
    n, T = S.shape
    for i in range(n):
        W = 0.0; Bg = 0.0; Wprev = 0.0; ip = 0.0; up = 0.0; hpr = 0.0
        for k in range(T):
            s = S[i, k] * 1.0
            hm = hs[k] * kt
            c = cnt[i, k]
            if c > 0.0:
                Bg = Bg + c
                Wk[i, k] = W; Wp[i, k] = Wprev; inc[i, k] = ip; unw[i, k] = up; hp[i, k] = hpr
            aV = np.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / isc
            wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
            Wb = W
            W = aV * W + bV * (wb * s)
            Bg = Bg * np.exp(-hm / tauD)
            ip = bV * wb * s; up = bV * s; hpr = hm; Wprev = Wb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--results", required=True)
    ap.add_argument("--protos", required=True)
    a = ap.parse_args()
    fit = json.load(open(a.fit)); fa = fit["args"]
    P = {**MV.DEFAULTS, **V5_DEFAULTS, **json.loads(fa["filters"]), **json.loads(fa.get("set", "{}")), **fit["pre"]}
    protos = a.protos.split(",")
    live = {}
    for l in open(a.results):
        x = json.loads(l)
        if x.get("ok") and x["path"] == "L5" and x["proto"] in protos and x["cond"] == "control":
            live[(x["pair"], x["proto"])] = {int(s): (x["end"]["necb_GB"][i], x["v7"]["nveto"][i]) for i, s in enumerate(x["syn"])}
    batch_v2.BASIS_DIR = C.BASIS["L5"]
    B = batch_v2.BatchV2(PATHS["L5"]["dirs"], protocols=sorted(protos), pairs={p for p, _ in live}, fast=False,
                         signals=("vdcc",))
    kt = 1e3 * float(batch_v2.TAU_IND_GB); thV = rho_v4.vamp(P)[1]
    Tv = float(P["veto_T"]); tauE1 = float(P["tau_E1"]); isc = float(P["i_scale"]); thE0 = float(P["theta_eCB"])
    um = live_v7.ue_map(fit, "L5")
    acc = {}
    for r in B.recs:
        if (r["pair"], r["proto"]) not in live:
            continue
        t = np.asarray(r["t"], np.float64); n, T = r["effcai"].shape
        hs = np.zeros(T); hs[:T - 1] = np.diff(t) / 1000.0 / batch_v2.TAU_IND_GB
        cnt = np.zeros((n, T)); np.add.at(cnt, (np.broadcast_to(np.arange(n)[:, None], r["arr"].shape), r["arr"]), 1.0)
        pre, post = map(int, r["pair"].split("-"))
        uE = np.ones(n) if um is None else np.array([um[(pre, post, int(s))] for s in r["syn"]])
        thE = uE * thE0 if thE0 > 0 else np.full(n, thV)
        Q = C7.veto_Q(r, cnt, Tv, tauE1, isc)
        S = np.ascontiguousarray(r["vdcc"], dtype=np.float32)
        Wk, Wp, inc, unw, hp = (np.zeros((n, T)) for _ in range(5))
        pools(S, hs, cnt, kt, isc, tauE1, float(P["tau_d_NMDA"]), Wk, Wp, inc, unw, hp)
        ta = np.asarray(r["prespikes"], float)[None, :] + np.asarray(MV.edge_params(r["syn"])["delay"])[:, None]
        post_sp = np.sort(np.asarray(r["postspikes"], float).ravel())
        A = acc.setdefault(r["proto"], dict(gap=[], dpost=[], wr=[], cnt={v: [0, 0, 0] for v in ("V0", "V1", "V2")},
                                            live=[0, 0], syn_eq={v: 0 for v in ("V0", "V1", "V2")}, nsyn=0))
        L = live[(r["pair"], r["proto"])]
        for i, s in enumerate(r["syn"]):
            k = np.asarray(r["arr"][i]).ravel(); tai = ta[i].ravel()
            g = t[k] - tai; f = np.where(hp[i, k] > 0, g / np.maximum(hp[i, k], 1e-12), 0.0)
            j = np.searchsorted(post_sp, tai) if len(post_sp) else None
            if len(post_sp):
                near = np.stack([post_sp[np.clip(j - 1, 0, len(post_sp) - 1)], post_sp[np.clip(j, 0, len(post_sp) - 1)]])
                dp = near - tai[None, :]; A["dpost"] += list(dp[np.argmin(np.abs(dp), axis=0), np.arange(len(tai))])
            A["gap"] += list(g); A["wr"] += list(Wk[i, k] / thE[i])
            V = {"V0": (Wk[i, k], Q[i, k]), "V1": (Wk[i, k] - f * inc[i, k], Q[i, k] + f * unw[i, k]),
                 "V2": (Wp[i, k], Q[i, k] + unw[i, k])}
            ln, lv = L.get(int(s), (np.nan, np.nan)); A["live"][0] += ln; A["live"][1] += lv; A["nsyn"] += 1
            for v, (w, q) in V.items():
                tr = w > thE[i]; st = int(np.sum(tr & ~(q > thE[i]))); ve = int(np.sum(tr & (q > thE[i])))
                A["cnt"][v][0] += int(tr.sum()); A["cnt"][v][1] += st; A["cnt"][v][2] += ve
                A["syn_eq"][v] += int(st == ln and ve == lv)
    for p, A in acc.items():
        q = lambda x: np.round(np.quantile(np.asarray(x), [0.0, 0.1, 0.5, 0.9, 1.0]), 3)
        wr = np.asarray(A["wr"])
        print(f"\n== {p}: {A['nsyn']} synapses, {len(A['gap'])} arrivals")
        print(f"  gap t_k - t_a (ms) q0/10/50/90/100 {q(A['gap'])}")
        if A["dpost"]:
            print(f"  nearest post spike - t_a (ms) {q(A['dpost'])}")
        print(f"  W(t_k)/thE {q(wr)}; arrivals with 0.5 < W/thE < 2: {int(np.sum((wr > 0.5) & (wr < 2)))}")
        print(f"  live: steps {A['live'][0]:.0f}, vetoes {A['live'][1]:.0f}, triggered {sum(A['live']):.0f}")
        for v, (tr, st, ve) in A["cnt"].items():
            print(f"  {v}: triggered {tr}, steps {st}, vetoes {ve}; synapses with (steps, vetoes) = live: "
                  f"{A['syn_eq'][v]}/{A['nsyn']}")


if __name__ == "__main__":
    main()
