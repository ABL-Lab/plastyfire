"""v8 fitter (S4WIDE_DIAG.md A/B): copy of fit_v7.py (not modified) with a final-point selection that replaces phase B.
  Phase B (the --resume call at the binary readout) used to be a near no-op: the collapsed phase-A population stops the DE after
  1-20 generations on a piecewise-constant objective. New default --phaseb rescore (only acts with --resume and --maxiter > 0):
    1. load the phase-A ckpt population (x, f smoothed energies), drop near-duplicate members, keep the top --phaseb-k (30) by f;
    2. score them with the binary readout (the DE objective incl. hinge); order by (binary, smoothed f), so plateau ties go to the
       smoother point;
    3. Nelder-Mead polish (scipy, in the normalised box, bounds respected, <= --polish-nfev evals each) from the best
       --polish-k (3) distinct points, on the same binary objective;
    4. the best of everything scored is the final point (json: de_fun = its binary objective, nfev = all evals).
  --phaseb de keeps the v7 behaviour (DE resume). Both modes add phaseA_{generations,nfev,minutes,de_fun} to the json, read from
  the ckpt (nit) and its sibling <ckpt minus _ckpt.npz>.json; the json also gets phaseB = {mode, n_scored, scored, polish}.
v7 fitter (SPEC_FITPROC.md, fit-procedure fixes 2 and 4, user-approved 2026-10-02): copy of fit_v6.py (not modified)
plus three flags. All off by default, so the code path and the results are those of fit_v6 (REPRO bit-identical).
  --init-admissible R  (default 0 = off): build the initial population as before (seeded: uniform + seeds in members 1..n;
      unseeded: a Latin hypercube from rng(--seed)), then replace every non-seed member that fails fit_v3.rules_row by an
      admissible candidate from up to R rounds of --init-pool x NP Latin-hypercube points. Seeds are kept as they are.
  --strategy S  (default best1bin): scipy differential_evolution strategy.
  --hinge JSON  (env HINGE, default off): {"NAME/target|cond": [k, lam], ...}; for those targets the DE objective uses
      lam * max(|z| - k, 0)^2 in place of w z^2 (NaN prediction: 1e3). The csv chi2 stays unweighted z^2.
Run through fit_launch.py (env FITTER=fit_v7) so a gpu_v1x kernel's `import fit_v6` loads this module.

v6 fitter (EXP_PARAM.md): copy of fit_v5.py (not modified) on gpu_v6_rho, which adds per-synapse threshold scales so
that the v5c thresholds can be written in the synapse's own measured Ca quantities (cexp csv of the per-synapse
Cpre / Cpre_APV / Cpre_Mg0 / Cpost measurement, columns pathway,pre_gid,post_gid,syn_id,loc,dist_um,cpre,cpre_apv,
cpre_mg0,cpost (+ *_cai, vdcc_q, nmda_q)).
  --cexp-dir DIR (env CEXP_DIR): <DIR>/<pathway>.csv joined onto every basis synapse by (pre_gid, post_gid, syn_id);
      pathway = CEXP_MAP[model] (default l5 -> l5l5, l23 -> l23l5, l23l23 -> l23l23); a missing synapse is an error.
      DIR = STANDIN: a clearly labelled stand-in table (cpre = npz c_pre, cpost = npz c_post, cpre_apv = 0.5 c_pre,
      cpre_mg0 = 2 c_pre, *_cai = the same numbers, loc STANDIN) is written to /scratch/dhuruva/split1/cexp_standin
      and joined like the real one (code-path smoke test only; its numbers mean nothing).
  --scale-V / --scale-E (env SCALE_V / SCALE_E): json {column: coef} -> u_i = sum coef * column_i; the kernel uses
      theta_V,i = theta_V u_V,i and theta_eCB,i = theta_eCB u_E,i, so the fitted theta_V / theta_eCB become the
      dimensionless k_V / k_E. Column "one" = 1. With 48.3 x a *_cai column, u is that event's pool level
      (C-E transient calibration, ANCHOR_V5C §2).
  --theta-pre (env THETA_PRE): json {"d": {col: coef}, "p": {col: coef}} replaces c_pre in theta_d / theta_p.
  Seeds: with scales active, the seed's theta_V / theta_eCB are divided by the median u (a00 / a10 by the median
      Xd / c_pre, Xp / c_pre) so the seed starts at the v5c thresholds of a median synapse.
  With CEXP real, the log prints the relabelling of the seed's thresholds (EXP_PARAM §1a): per pathway, quantiles of
      theta / (48.3 C_cai) for each measured event, and writes the per-synapse table to /scratch/dhuruva/split1/cexp_relabel.
  No scales and no theta-pre: the rule is v5 exactly; --maxiter 0 on a v5/v6 json asserts the repro (1e-6).
Everything below is fit_v5.py.

v5 fitter (V5_DESIGN.md): fit_v4n.py (N pathway models, --extra) with the v5 shared-pool rule of gpu_v5_rho.
fit_v4n.py / fit_v4.py are not modified.

v5 = Chindemi 2022 + one own-spine VDCC-Ca pool V (tau_E1, i_scale) with one threshold theta_V that both gates
potentiation and lets an own pre spike drive eCB-LTD (d -> d_min). No NO arm, rho_gamma 1. Set by the filters
(run_fit_v5.sh): v5_mode 1, vamp_mode 1, rho_gamma 1, A_eCB, dpre_min, tau_E1, i_scale. v5_mode 2 = v5b ((1 - b)-weighted
eCB trigger pool, threshold theta_V); + --free-filters theta_eCB = v5c (separate eCB threshold, log box 0.05-50);
v5_mode 3 = v5-E2 (v5b with only VDCC excursions that start in an unprimed pool, V <= theta_V, feeding the trigger).
DE vector: a00, a01, a10, a11 + --free-filters (theta_V) + gamma_d, gamma_p (--fit-gamma). The v4 A_mglu / A_NO slots
are not in the vector: they are held at 0 exactly (log slot -inf, fit_v2.unpack gives 10**-inf = 0.0), so the v4
eCB / NO chain is off in every v4 evaluation made here.

--check-v4: before the DE, at the first seed (a's, gammas, theta_V of the seed; rho_gamma 1, A_mglu = A_NO = 0), for
every pathway model assert (1e-9 in chi2)
  M1: v5 kernel (A_eCB 0)               == v4 vamp kernel (vamp_mode 1, gpu_v4_rho)
  M0: v5 kernel (A_eCB 0, theta_V 0)    == v4 rates kernel (vamp_mode 0) == v4 vamp kernel at theta_V 0
and print the v4 M0 chi2 next to results/v4_LM0chk.json (ladder check 22135018) when that file exists.
With --maxiter 0 the output is the seed scored with the full v5 rule (the filters).

--extra NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]] (repeatable) and --joint as in fit_v4n.py (one uniform parameter set;
pathway differences only via each synapse's own traces and basis). --drop-targets "target|cond" (all models) or
"NAME/target|cond". Outputs: <save>.csv (l5), <save>_<NAME>.csv, json chi2_<NAME>, n_targets_<NAME>, chi2_total.
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import fit_v2 as F                    # noqa: E402
import fit_v3 as F3                   # noqa: E402  (rules_row)
import batch_v2                       # noqa: E402
from batch_v2 import BatchV2          # noqa: E402
from targets import load_targets      # noqa: E402
import model_v2 as MV                 # noqa: E402
import rho_v4                         # noqa: E402
import gpu_v5_rho                     # noqa: E402
import gpu_v6_rho                     # noqa: E402

F.FILTER_BOX.update(rho_gamma=(0.0, 1.0), tau_fast=(5.0, 300.0), rho_sigma=(0.005, 0.3),   # runtime only
                    gamma_d=(20.0, 250.0), gamma_p=(100.0, 600.0),    # GAMMA_LIT.md, as fit_v4n
                    theta_V=(0.0, 50.0),                              # shared pool threshold (flagged fitted)
                    theta_eCB=(0.05, 50.0))                           # v5c separate eCB threshold on W (log box)
F.LOG_FILTERS |= {"tau_fast", "rho_sigma", "theta_eCB"}
for k, v in {**rho_v4.V4_DEFAULTS, **rho_v4.RATE_DEFAULTS, **rho_v4.VAMP_DEFAULTS, **gpu_v5_rho.V5_DEFAULTS}.items():
    MV.DEFAULTS.setdefault(k, v)

GEOM = os.environ.get("GEOM_L23") or os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")   # GEOM_L23: alternative geometry csv
L5_BASIS = batch_v2.BASIS_DIR
L23_BASIS = os.environ.get("L23_BASIS_DIR") or os.path.join(ROOT, "basis_results_edges_ebner_l23l5_delta_rs")   # L23_BASIS_DIR: other emodel (run_new_emodel_chain.sh)
LM0CHK = os.path.join(HERE, "results", "v4_LM0chk.json")
_M = []            # [(model, name)]; _M[0] is L5 (rules use its peaks)
_CFG = None
_POP = {}


def expand(x):
    """DE vector (no A_mglu / A_NO) -> fit_v2 vector with both log slots at -inf (amplitudes exactly 0)."""
    return np.concatenate([x[:4], [-np.inf, -np.inf], x[4:]])


def unpack5(x):
    a, P = F.unpack(expand(np.asarray(x, float)), _CFG)
    return a, P


def build(dirs, groups, conds, fil, pairs, basis, fast, loc=None, vgate=False, rates=True, drop=()):
    batch_v2.BASIS_DIR = basis
    T = {k: v for k, v in load_targets(tuple(groups.split(","))).items() if k[1] in conds and f"{k[0]}|{k[1]}" not in drop}
    protos = sorted({k[0].split("@")[0] for k in T})
    sig = ("vdcc",) if fil.get("pre_drive") else ("shaft_cai",)
    B = BatchV2(dirs.split(","), protocols=protos, pairs=pairs, fast=False, signals=sig)
    have = {r["proto"] for r in B.recs}
    T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    peak = np.concatenate([r["peak"] for r in B.recs])
    B._stack()
    G = gpu_v6_rho.GPUModelV6(B, T, fil, fast=fast, loc=loc, vgate=vgate, rates=rates)
    G._keys = [(r["pair"], int(s)) for r in B.recs for s in r["syn"]]          # kernel order = G.cpre order
    assert len(G._keys) == G.n_ctrl == len(G.cpre)
    for r in B.recs:
        for k in ("effcai", "vdcc", "shaft_cai", "cacr", "t", "cev", "cev_lo", "cev_hi", "vev", "_c4", "_v3", "arr"):
            r.pop(k, None)
    B._E = B._H = None
    import gc; gc.collect()
    print(f"{groups}: {len(T)} targets, {len(B.recs)} records (v5 kernel {G.v5})", flush=True)
    return G, T, peak


CEXP_MAP = {"l5": "L5L5", "l23": "L23L5", "l23l23": "L23L23", **json.loads(os.environ.get("CEXP_MAP") or "{}")}   # A4 file names
STANDIN_DIR = "/scratch/dhuruva/split1/cexp_standin"
RELABEL_DIR = "/scratch/dhuruva/split1/cexp_relabel"
POOL_PER_UM = 48.3        # pool units per uM peak free spine Ca of one event (C-E transient calibration, ANCHOR_V5C §2)
EVENTS = ("cpre", "cpre_apv", "cpre_mg0", "cpost")


def write_standin(G, name):
    """STAND-IN cexp table from the basis npz c_pre / c_post (labelled loc STANDIN). Smoke test of the join only."""
    import pandas as pd
    seen, rows = set(), []
    for (pair, s), cp, cq in zip(G._keys, G.cpre, G.cpost):
        if (pair, s) in seen:
            continue
        seen.add((pair, s)); pre, post = pair.split("-")
        rows.append(dict(pathway=CEXP_MAP[name], pre_gid=int(pre), post_gid=int(post), syn_id=s, loc="STANDIN",
                         dist_um=np.nan, cpre=cp, cpre_apv=0.5 * cp, cpre_mg0=2.0 * cp, cpost=cq, cpre_cai=cp,
                         cpre_apv_cai=0.5 * cp, cpre_mg0_cai=2.0 * cp, cpost_cai=cq, vdcc_q=np.nan, nmda_q=np.nan))
    os.makedirs(STANDIN_DIR, exist_ok=True)
    f = os.path.join(STANDIN_DIR, f"{CEXP_MAP[name]}.csv")
    pd.DataFrame(rows).to_csv(f, index=False)
    print(f"STAND-IN cexp (not a measurement): {len(rows)} synapses -> {f}", flush=True)


def join_cexp(G, name, cdir):
    """cexp rows in kernel order (DataFrame, len n_ctrl). Every basis synapse must be in the table."""
    import pandas as pd
    f = os.path.join(cdir, f"{CEXP_MAP[name]}.csv")
    df = pd.read_csv(f)
    df["key"] = df.pre_gid.astype(str) + "-" + df.post_gid.astype(str) + ":" + df.syn_id.astype(np.int64).astype(str)
    if df.key.duplicated().any():
        raise SystemExit(f"{f}: {int(df.key.duplicated().sum())} duplicated synapse keys")
    keys = [f"{p}:{s}" for p, s in G._keys]
    idx = df.set_index("key").reindex(keys)
    miss = idx.cpre.isna()
    if miss.any():
        raise SystemExit(f"{f}: {int(miss.sum())} of {len(keys)} basis synapse rows not in cexp "
                         f"(first {[k for k, m in zip(keys, miss) if m][:3]})")
    idx = idx.reset_index()
    for c, ref in (("cpre", G.cpre), ("cpost", G.cpost)):
        rel = np.abs(idx[c].to_numpy(float) - ref) / np.maximum(np.abs(ref), 1e-12)
        print(f"cexp {name} ({f}): {c} vs npz c_{c[1:]}: rel diff q50 {np.nanmedian(rel):.2e} max {np.nanmax(rel):.2e}, "
              f"NaN in cexp {int(np.isnan(rel).sum())} of {len(rel)} rows", flush=True)
    return idx


def lincomb(df, spec, N, base=None):
    """sum coef * column. Column "one" = 1; "a/b" = df[a] / df[b] (cexp ratio), times base when given (theta-pre:
    base = npz c_pre, so only the cexp ratio enters and the Chindemi c_pre stays the cache value)."""
    u = np.zeros(N)
    for col, coef in spec.items():
        if col == "one":
            v = np.ones(N)
        elif "/" in col:
            a, b = col.split("/")
            v = df[a].to_numpy(np.float64) / df[b].to_numpy(np.float64)
            if base is not None:
                v = v * base
        else:
            v = df[col].to_numpy(np.float64)
        u = u + float(coef) * v
    return u


def relabel(G, name, df, thV, thE, tag):
    """EXP_PARAM §1a: the seed's absolute thresholds as multiples of each measured event's pool level, per unique
    synapse. Two bridges: cai (48.3 x peak spine Ca, ANCHOR_V5C §2) and q (own VDCC charge of the event / i_scale,
    the pool's own input: 1e5 x vdcc_q[nA ms], peak pool level with no decay)."""
    df = df.drop_duplicates(["pre_gid", "post_gid", "syn_id"]).reset_index(drop=True)
    out = df[[c for c in ("pathway", "pre_gid", "post_gid", "syn_id", "loc", "dist_um") if c in df]].copy()
    q = lambda x: " ".join(f"{v:.3g}" for v in np.nanquantile(np.asarray(x, float), [0.1, 0.25, 0.5, 0.75, 0.9]))
    fr = lambda x: float(np.mean(np.asarray(x, float)[np.isfinite(x)] < 1)) if np.isfinite(x).any() else np.nan
    QSFX = dict(cpre="", cpre_apv="_apv", cpre_mg0="_mg0", cpost="_post")
    for ev in EVENTS:
        for br, c, conv in (("cai", f"{ev}_cai", POOL_PER_UM), ("q", f"vdcc_q{QSFX[ev]}", 1e5)):
            if c not in df:
                continue
            P = conv * df[c].to_numpy(float)
            with np.errstate(divide="ignore", invalid="ignore"):
                rv = thV / P; re = thE / P
            out[f"P_{br}_{ev}"] = P; out[f"thV_over_{br}_{ev}"] = rv; out[f"thE_over_{br}_{ev}"] = re
            print(f"relabel {name} {ev} [{br}: {conv:g} x {c}] n {len(P)} (NaN {int(np.isnan(P).sum())}): pool of one "
                  f"event q10/25/50/75/90 {q(P)} | theta_V {thV:.4g} / P: {q(rv)} (frac<1 {fr(rv):.2f}) | theta_eCB "
                  f"{thE:.4g} / P: {q(re)} (frac<1 {fr(re):.2f})", flush=True)
    cp, ca, cm = (df[k].to_numpy(float) for k in ("cpre", "cpre_apv", "cpre_mg0"))
    with np.errstate(divide="ignore", invalid="ignore"):
        print(f"relabel {name}: effcai NMDA share (Cpre - Cpre_APV)/Cpre {q((cp - ca) / cp)}; Mg-block depth "
              f"(Cpre_Mg0 - Cpre)/Cpre {q((cm - cp) / cp)}; Cpost/Cpre {q(df.cpost.to_numpy(float) / cp)}", flush=True)
    os.makedirs(RELABEL_DIR, exist_ok=True)
    f = os.path.join(RELABEL_DIR, f"{tag}_{name}.csv"); out.to_csv(f, index=False)
    print(f"relabel {name}: {len(out)} unique synapses -> {f}", flush=True)


def objective_gpu(X):
    xs = X.T if X.ndim == 2 else X[None]
    out = np.zeros(len(xs))
    AP = [unpack5(x) for x in xs]
    Ps_all = [{**_CFG["filters"], **P} for _, P in AP]
    G0 = _M[0][0]
    td, tp = G0.thetas_all([a for a, _ in AP], Ps_all)
    idx = []
    for i, (a, P) in enumerate(AP):
        pen = F3.rules_row(td[i], tp[i], _CFG["peak"])
        if pen:
            out[i] = pen
        else:
            idx.append(i)
    if idx:
        t0 = time.time()
        Ps = [Ps_all[i] for i in idx]; A = [AP[i][0] for i in idx]
        c = G0.evaluate_th(td[idx], tp[idx], Ps)
        for G, _ in _M[1:]:
            c = c + G.evaluate(A, Ps)
        out[idx] = np.where(np.isfinite(c), c, 1e5)
        _POP["sec"] = time.time() - t0; _POP["n"] = len(idx)
    return out if X.ndim == 2 else out[0]


def rules_pen(X, chunk=512):
    """v7: fit_v3.rules_row penalty of each row of X (n, D) on the L5 model, in chunks (0 = admissible)."""
    out = np.empty(len(X))
    for s in range(0, len(X), chunk):
        AP = [unpack5(x) for x in X[s:s + chunk]]
        td, tp = _M[0][0].thetas_all([a for a, _ in AP], [{**_CFG["filters"], **P} for _, P in AP])
        out[s:s + len(AP)] = [F3.rules_row(td[i], tp[i], _CFG["peak"]) for i in range(len(AP))]
    return out


WEIGHTS = {}       # env WEIGHTS / --weights: {"ltd": w, "ltp": w, "target|cond": w, "NAME/target|cond": w}; empty = unweighted
HINGE = {}         # v7 --hinge: {"target|cond" or "NAME/target|cond": [k, lam]}; empty = off


def target_weight(name, k, mean):
    """z^2 weight of one target: ltd (exp mean < 1.0) / ltp class weight, overridden by an explicit key."""
    w = WEIGHTS.get("ltd", 1.0) if mean < 1.0 else WEIGHTS.get("ltp", 1.0)
    for key in (f"{k[0]}|{k[1]}", f"{name}/{k[0]}|{k[1]}"):
        w = WEIGHTS.get(key, w)
    return float(w)


def hinge_for(name, T):
    """v7: {target key: (k, lam)} of this model's targets named in HINGE (plain or NAME/ prefixed key)."""
    H = {}
    for k in T:
        for key in (f"{k[0]}|{k[1]}", f"{name}/{k[0]}|{k[1]}"):
            if key in HINGE:
                H[k] = (float(HINGE[key][0]), float(HINGE[key][1]))
    return H


def install_weights(G, T, name):
    """Instance-level chi2 wrapper: weighted sum of z^2 for the DE objective (class chi2 untouched; only called when
    WEIGHTS is set, so the default path is bit-identical)."""
    orig = type(G).chi2; W = {k: target_weight(name, k, T[k][0]) for k in T}
    H = hinge_for(name, T)

    def chi2w(rho, d, return_pred=False, Ps=None):
        _, preds = orig(G, rho, d, True, Ps)
        out = np.zeros(rho.shape[0])
        for k, v in preds.items():
            z = (v - T[k][0]) / T[k][1]
            if k in H:                                  # v7 hinge: replaces z^2 and its weight for this target
                kk, lam = H[k]
                out += lam * np.where(np.isfinite(z), np.maximum(np.abs(z) - kk, 0.0) ** 2, 1e3)
            else:
                out += W[k] * np.where(np.isfinite(z), z ** 2, 0.0)
        return (out, preds) if return_pred else out
    G.chi2 = chi2w


def table(G, T, a, P, name="l5"):
    """Unweighted chi2 (class chi2, also when a weighting wrapper is installed) + per-target table."""
    import pandas as pd
    c, preds = type(G).chi2(G, *G.rho_dpre([a], [P]), return_pred=True, Ps=[P])
    df = pd.DataFrame([dict(target=k[0], condition=k[1], target_mean=T[k][0], target_sem=T[k][1],
                            pred=float(v[0]), z=(float(v[0]) - T[k][0]) / T[k][1]) for k, v in preds.items()])
    if WEIGHTS:
        df["ltd"] = df.target_mean < 1.0
        df["w"] = [target_weight(name, (t, c_), m) for t, c_, m in zip(df.target, df.condition, df.target_mean)]
    return float(c[0]), df


def check_v4(seed, fil):
    """M1 / M0 limits of the v5 kernel against the v4 kernels at the seed (see module doc). Raises on a mismatch."""
    import gpu_v4_rho
    f, fj, x = seed
    a, P = unpack5(x); P = {**fil, **P}
    assert P["A_mglu"] == 0.0 and P["A_NO"] == 0.0, (P["A_mglu"], P["A_NO"])
    P1 = {**P, "A_eCB": 0.0}                      # M1: C1 gate on the pool, no eCB, no NO
    P0 = {**P1, "theta_V": 0.0}                   # M0: pure Chindemi (fitted rates)
    tot = {"m1_v4": 0.0, "m1_v5": 0.0, "m0_v4": 0.0, "m0_v5": 0.0, "v5": 0.0}
    for G, name in _M:
        assert G.v5 and G.vamp == 1 and G.rates, (name, G.v5, G.vamp, G.rates)
        kv = G._kern
        kr = gpu_v4_rho._make_kernel_fast(G.t_mode, G.no_mode, G.use_s, False, False, True, 0)   # v4 rates, vamp 0

        def ev(v5, Pq, kern=None):
            G.v5 = v5
            if kern is not None:
                G._kern = kern
            try:
                return float(G.evaluate([a], [Pq])[0])
            finally:
                G._kern = kv; G.v5 = True
        c = dict(m1_v4=ev(False, P1), m1_v5=ev(True, P1), m0_v4=ev(False, P0, kr), m0_v4g=ev(False, P0),
                 m0_v5=ev(True, P0), v5=ev(True, P))
        for k in tot:
            tot[k] += c[k]
        print(f"check-v4 {name}: M1 v4 {c['m1_v4']:.12f}  v5 {c['m1_v5']:.12f}  diff {abs(c['m1_v4'] - c['m1_v5']):.1e} | "
              f"M0 v4 {c['m0_v4']:.12f}  v4-vamp(theta_V 0) {c['m0_v4g']:.12f}  v5 {c['m0_v5']:.12f}  "
              f"diff {abs(c['m0_v4'] - c['m0_v5']):.1e} | v5 full rule at seed {c['v5']:.6f}", flush=True)
        assert abs(c["m1_v4"] - c["m1_v5"]) < 1e-9, (name, "M1", c)
        assert abs(c["m0_v4"] - c["m0_v5"]) < 1e-9 and abs(c["m0_v4"] - c["m0_v4g"]) < 1e-9, (name, "M0", c)
    print(f"check-v4 total over {len(_M)} models: M1 v4 {tot['m1_v4']:.9f} v5 {tot['m1_v5']:.9f}; "
          f"M0 v4 {tot['m0_v4']:.9f} v5 {tot['m0_v5']:.9f}; v5 full {tot['v5']:.6f}  (seed {os.path.basename(f)}: "
          f"json chi2_total {fj.get('chi2_total', fj['chi2'] + fj.get('chi2_l23', 0.0)):.6f})", flush=True)
    if os.path.exists(LM0CHK):
        lj = json.load(open(LM0CHK))
        ref = lj.get("chi2", np.nan) + lj.get("chi2_l23", 0.0)
        print(f"check-v4 vs ladder {os.path.basename(LM0CHK)}: json total {ref:.9f}  v4 M0 here {tot['m0_v4']:.9f}  "
              f"diff {abs(ref - tot['m0_v4']):.1e} (equal only if the ladder used these a's / gammas)", flush=True)
    else:
        print(f"check-v4: {LM0CHK} not there yet (ladder check pending), in-job v4 kernels used", flush=True)
    print("CHECK-V4 OK (1e-9)", flush=True)
    return tot


def rescore_polish(ck, lo, hi, args):
    """v8 phase B: binary rescoring of the phase-A population + Nelder-Mead polish. Returns (x, binary objective, nfev, info)."""
    from scipy.optimize import minimize
    X = np.asarray(ck["x"], float); f = np.asarray(ck["f"], float); span = hi - lo
    U = (X - lo) / span
    order = np.argsort(f, kind="stable"); keep = []
    for i in order:                                  # distinct members (normalised L_inf > 1e-4), best phase-A energy first
        if all(np.max(np.abs(U[i] - U[j])) > 1e-4 for j in keep):
            keep.append(i)
        if len(keep) >= args.phaseb_k:
            break
    keep = np.array(keep)
    t1 = time.time()
    bin_f = np.asarray(objective_gpu(X[keep].T), float)
    nfev = len(keep)
    srt = np.lexsort((f[keep], bin_f))               # binary first, smoothed energy breaks plateau ties
    keep, bin_f = keep[srt], bin_f[srt]
    print(f"phase B rescore: {len(keep)} distinct of {len(X)} members scored at the binary readout in {time.time() - t1:.1f}s; "
          f"binary best {bin_f[0]:.4f} median {np.median(bin_f):.4f} worst {bin_f[-1]:.4f}; phase-A f of the binary best "
          f"{f[keep[0]]:.4f} (phase-A min {f.min():.4f})", flush=True)
    best_x, best_f = X[keep[0]].copy(), float(bin_f[0])
    info = dict(mode="rescore", n_members=len(X), n_scored=len(keep), scored_binary=[float(v) for v in bin_f],
                scored_phaseA=[float(f[i]) for i in keep], rescored_best=best_f, polish=[])
    starts = []                                      # best distinct (normalised L_inf > 0.02) points
    for i in keep:
        if all(np.max(np.abs(U[i] - U[j])) > 0.02 for j in starts):
            starts.append(i)
        if len(starts) >= args.polish_k:
            break
    D = len(lo)
    for i in (starts if args.polish_nfev > 0 else []):
        u0 = U[i]; n = [0]; seen = [float("inf"), u0]; t2 = time.time()

        def fo(u):
            n[0] += 1; v = float(objective_gpu(lo + np.clip(u, 0.0, 1.0) * span))
            if v < seen[0]:
                seen[0], seen[1] = v, np.clip(u, 0.0, 1.0)
            return v
        sim = np.vstack([u0] + [np.clip(u0 + np.where(np.arange(D) == k, -0.03 if u0[k] > 0.5 else 0.03, 0.0), 0, 1)
                                for k in range(D)])
        r = minimize(fo, u0, method="Nelder-Mead", bounds=[(0.0, 1.0)] * D,
                     options=dict(maxfev=args.polish_nfev, initial_simplex=sim, xatol=1e-5, fatol=1e-9))
        nfev += n[0]; xp = lo + seen[1] * span
        info["polish"].append(dict(start_binary=float(bin_f[list(keep).index(i)]), end_binary=seen[0], nfev=n[0],
                                   sec=time.time() - t2))
        print(f"polish from member {i}: binary {info['polish'][-1]['start_binary']:.4f} -> {seen[0]:.4f} in {n[0]} evals "
              f"({time.time() - t2:.0f}s)", flush=True)
        if seen[0] < best_f:
            best_x, best_f = xp, seen[0]
    info["final_binary"] = best_f
    print(f"phase B final: binary objective {best_f:.4f} (rescored best {info['rescored_best']:.4f}, nfev {nfev})", flush=True)
    return best_x, best_f, nfev, info


def pd_read(f):
    import pandas as pd
    return pd.read_csv(f)


def main():
    global _CFG
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", required=True)
    ap.add_argument("--groups", default="paired_l5,sjostrom07")
    ap.add_argument("--filters", default="{}")
    ap.add_argument("--set", default="{}", help="options fixed for every candidate (merged into filters)")
    ap.add_argument("--free-filters", default="theta_V")
    ap.add_argument("--conditions", default="control,mglu_block,post_nmdar,nmdar_block,no_block")
    ap.add_argument("--a-lo", type=float, default=0.0); ap.add_argument("--a-hi", type=float, default=5.0)
    ap.add_argument("--x0", default=None)
    ap.add_argument("--maxiter", type=int, default=150)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", default=None)
    ap.add_argument("--save", required=True)
    ap.add_argument("--seed-fits", default=None)
    ap.add_argument("--seed-set", default="{}")
    ap.add_argument("--resume", default=None)
    ap.add_argument("--joint", action="store_true")
    ap.add_argument("--l23-dirs", default=os.path.join(V2, "extracted", "ebner_l23l5_delta-prefire-vseg-rs"))
    ap.add_argument("--check-v4", action="store_true")
    ap.add_argument("--drop-targets", default="", help="target|condition,... left out of the objective (validation only)")
    ap.add_argument("--fit-gamma", action="store_true", help="free gamma_d, gamma_p (rho rates, uniform)")
    ap.add_argument("--weights", default=os.environ.get("WEIGHTS", ""),
                    help='json {"ltd": w, "ltp": w, "target|cond": w}: z^2 weights in the DE objective only (env WEIGHTS)')
    ap.add_argument("--extra", action="append", default=[], help="NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]], repeatable")
    ap.add_argument("--cexp-dir", default=os.environ.get("CEXP_DIR", ""), help="cexp csv dir or STANDIN (env CEXP_DIR)")
    ap.add_argument("--scale-V", default=os.environ.get("SCALE_V", ""), help='json {col: coef}: theta_V,i = theta_V u_i')
    ap.add_argument("--scale-E", default=os.environ.get("SCALE_E", ""), help='json {col: coef}: theta_eCB,i = theta_eCB u_i')
    ap.add_argument("--theta-pre", default=os.environ.get("THETA_PRE", ""), help='json {"d": {col: coef}, "p": {...}}')
    ap.add_argument("--box", default=os.environ.get("BOX", ""), help='json {filter: [lo, hi]} overrides FILTER_BOX')
    ap.add_argument("--init-admissible", type=int, default=0,
                    help="v7: max LHS rounds to replace inadmissible non-seed initial members (0 = off, fit_v6 path)")
    ap.add_argument("--init-pool", type=int, default=20, help="v7: candidates per round = init_pool x NP")
    ap.add_argument("--strategy", default="best1bin", help="v7: differential_evolution strategy")
    ap.add_argument("--hinge", default=os.environ.get("HINGE", ""),
                    help='v7: json {"NAME/target|cond": [k, lam]}: lam max(|z| - k, 0)^2 replaces w z^2 in the DE objective')
    ap.add_argument("--phaseb", default="rescore", choices=["rescore", "de"],
                    help="v8: with --resume, rescore the phase-A population at the binary readout + Nelder-Mead polish "
                         "(default) or de = the v7 DE resume")
    ap.add_argument("--phaseb-k", type=int, default=30, help="v8: members (top by phase-A energy, distinct) to rescore")
    ap.add_argument("--polish-k", type=int, default=3, help="v8: best distinct points that get a local polish (0 = none)")
    ap.add_argument("--polish-nfev", type=int, default=200, help="v8: Nelder-Mead evaluations per polished point")
    args = ap.parse_args()
    if args.hinge:
        HINGE.update({k: tuple(v) for k, v in json.loads(args.hinge).items()})
    if args.box:
        F.FILTER_BOX.update({k: tuple(v) for k, v in json.loads(args.box).items()})
    sV, sE = (json.loads(s) if s else {} for s in (args.scale_V, args.scale_E))
    tpre = json.loads(args.theta_pre) if args.theta_pre else {}
    if (sV or sE or tpre) and not args.cexp_dir and any(c != "one" for s in (sV, sE, *tpre.values()) for c in s):
        raise SystemExit("--scale-V / --scale-E / --theta-pre with measured columns need --cexp-dir")

    conds = set(args.conditions.split(","))
    fil = {**json.loads(args.filters), **json.loads(args.set)}
    free = [k for k in args.free_filters.split(",") if k]
    if args.fit_gamma:
        free += [k for k in ("gamma_d", "gamma_p") if k not in free]
    specs = []                                   # (name, groups, dirs, basis, geom, pairs_file)
    if args.joint:
        specs.append(("l23", "paired_l23l5", args.l23_dirs, L23_BASIS, GEOM, ""))
    for s in args.extra:
        f = s.split(":")
        if len(f) < 4 or len(f) > 6:
            raise SystemExit(f"--extra {s!r}: need NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]]")
        f += [""] * (6 - len(f))
        if not os.path.isabs(f[3]) and not os.path.isdir(f[3]):
            f[3] = os.path.join(ROOT, f[3])
        specs.append(tuple(f))
    names = ["l5"] + [s[0] for s in specs]
    if len(set(names)) != len(names):
        raise SystemExit(f"model names must be unique: {names}")
    drops = [t for t in args.drop_targets.split(",") if t]

    def drop_for(name):
        return set(t.split("/", 1)[1] if "/" in t else t for t in drops if "/" not in t or t.split("/", 1)[0] == name)

    G5, T5, peak = build(args.dirs, args.groups, conds, fil, set(args.pairs.split(",")) if args.pairs else None,
                         L5_BASIS, False, drop=drop_for("l5"))
    _M.append((G5, "l5"))
    TM = {"l5": T5}
    for name, groups, dirs, basis, geom, pf in specs:
        loc = pairs = None
        if geom:
            import pandas as pd
            g = pd.read_csv(geom); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
            lc = [c for c in g.columns if c.endswith("_distal")]
            loc = {p: {c: bool(v) for c, v in zip(lc, row)} for p, row in zip(g.pair, g[lc].itertuples(index=False))}
            pairs = set(g.pair)
        if pf:
            pp = set(open(pf).read().strip().split(","))
            pairs = pp if pairs is None else pairs & pp
        G, T, _ = build(dirs, groups, conds, fil, pairs, basis, False, loc=loc, drop=drop_for(name))
        _M.append((G, name)); TM[name] = T
        print(f"model {name}: groups {groups}, basis {basis}, geom {geom or '-'}, {len(T)} targets", flush=True)
    if args.weights or HINGE:
        if args.weights:
            WEIGHTS.update(json.loads(args.weights))
        for G, name in _M:
            install_weights(G, TM[name], name)
        for G, name in _M:
            ws = [target_weight(name, k, TM[name][k][0]) for k in TM[name]]
            print(f"weights {name}: {WEIGHTS}; {sum(w != 1.0 for w in ws)} of {len(ws)} targets weighted", flush=True)
            if HINGE:
                print(f"hinge {name}: {hinge_for(name, TM[name])} (lam max(|z| - k, 0)^2 replaces w z^2)", flush=True)
        if HINGE and not any(hinge_for(n, TM[n]) for _, n in _M):
            raise SystemExit(f"--hinge {HINGE}: no key matches a fitted target")
    # ---- v6: per-synapse measured quantities and threshold scales ---------------------------------------------------
    cdir = args.cexp_dir
    if cdir == "STANDIN":
        for G, name in _M:
            write_standin(G, name)
        cdir = STANDIN_DIR
    CX = {name: (join_cexp(G, name, cdir) if cdir else None) for G, name in _M}
    med = dict(uV=1.0, uE=1.0, xd=1.0, xp=1.0)
    if sV or sE or tpre:
        # Undefined measurement (NaN, e.g. cpost / vdcc_q_post at post cells that cannot fire one AP): that synapse
        # gets the median-synapse value (u = pooled median of the finite u over all models; Xd/Xp = c_pre x pooled
        # median of X / c_pre), i.e. exactly the v5c absolute threshold there. Same rule for every pathway.
        raw = {}
        for G, name in _M:
            N = G.n_ctrl; df = CX[name] if CX[name] is not None else {}
            cpre = np.maximum(G.cpre, 1e-12)
            raw[name] = dict(uV=lincomb(df, sV, N) if sV else None, uE=lincomb(df, sE, N) if sE else None,
                             xd=lincomb(df, tpre["d"], N, G.cpre) / cpre if "d" in tpre else None,
                             xp=lincomb(df, tpre["p"], N, G.cpre) / cpre if "p" in tpre else None)
        for k in med:
            v = [raw[n][k] for _, n in _M if raw[n][k] is not None]
            if v:
                v = np.concatenate(v); med[k] = float(np.median(v[np.isfinite(v)]))
        for G, name in _M:
            R = raw[name]; nan = {}
            for k in R:
                if R[k] is not None:
                    bad = ~np.isfinite(R[k]); nan[k] = int(bad.sum()); R[k] = np.where(bad, med[k], R[k])
            G.set_scales(R["uV"], R["uE"], None if R["xd"] is None else R["xd"] * G.cpre,
                         None if R["xp"] is None else R["xp"] * G.cpre)
            print(f"v6 scales {name}: uV {'-' if R['uV'] is None else np.quantile(R['uV'], [0.1, 0.5, 0.9]).round(4)}  "
                  f"uE {'-' if R['uE'] is None else np.quantile(R['uE'], [0.1, 0.5, 0.9]).round(4)}  "
                  f"Xd/c_pre {'-' if R['xd'] is None else np.quantile(R['xd'], [0.1, 0.5, 0.9]).round(4)}  "
                  f"Xp/c_pre {'-' if R['xp'] is None else np.quantile(R['xp'], [0.1, 0.5, 0.9]).round(4)}  "
                  f"NaN -> median rows {nan} of {G.n_ctrl}", flush=True)
        print(f"v6 seed conversion (row medians over all models): {med}  scale_V {sV} scale_E {sE} theta_pre {tpre}",
              flush=True)
    _CFG = dict(model="v2", rho="full", filters=fil, free_filters=free, filter_steps=8, loc=None, peak=peak,
                tie={}, continuous=True)
    F._CFG = _CFG

    bounds = [(args.a_lo, args.a_hi)] * 4 + [(0.0, 1.0)] * len(free)
    lo, hi = np.array(bounds).T
    print(f"v6 free ({len(bounds)}): {F.NAMES_A + free}; fixed filters {fil}", flush=True)
    seeds = []
    for f in (args.seed_fits.split(",") if args.seed_fits else []):
        fj = json.load(open(f)); pre = {**json.loads(args.seed_set), **fj["pre"]}
        aa = dict(fj["a"])
        if not fj.get("v6", {}).get("scaled"):          # absolute (v5 / unscaled) seed -> scaled units of a median synapse
            if sV and "theta_V" in pre:
                pre["theta_V"] = pre["theta_V"] / med["uV"]
            if sE and pre.get("theta_eCB", 0.0) > 0:
                pre["theta_eCB"] = pre["theta_eCB"] / med["uE"]
            aa["a00"] = aa["a00"] / med["xd"]; aa["a10"] = aa["a10"] / med["xp"]
        if cdir and not seeds and not fj.get("v6", {}).get("scaled"):
            for G, name in _M:
                relabel(G, name, CX[name], float(fj["pre"].get("theta_V", 0.0)), float(fj["pre"].get("theta_eCB", 0.0)),
                        os.path.basename(args.save) + ("_STANDIN" if args.cexp_dir == "STANDIN" else ""))
        with np.errstate(divide="ignore"):
            xf = F.pack(aa, pre, _CFG)
        seeds.append((f, fj, np.clip(np.delete(xf, [4, 5]), lo, hi)))
    if args.check_v4:
        assert seeds, "--check-v4 needs --seed-fits"
        check_v4(seeds[0], fil)
    init = "latinhypercube"
    if args.x0 or seeds:
        rng = np.random.default_rng(args.seed)
        pop = lo + rng.random((args.popsize * len(bounds), len(bounds))) * (hi - lo)
        if args.x0:
            a0 = json.load(open(args.x0)); a0 = a0.get("params", a0); pop[0, :4] = [a0[k] for k in F.NAMES_A]
        for i, (f, fj, x) in enumerate(seeds):
            pop[1 + i] = x; print(f"seed member {1 + i}: {f} (chi2 {fj['chi2']:.2f})", flush=True)
        init = pop
    if args.init_admissible > 0 and not args.resume:
        # v7: admissible initial population (SPEC_FITPROC.md fix 2); seed / x0 rows are kept as they are
        from scipy.stats import qmc
        rng = np.random.default_rng(args.seed); D = len(bounds); NP = args.popsize * D
        if isinstance(init, str):
            pop = lo + qmc.LatinHypercube(d=D, seed=rng).random(NP) * (hi - lo); fixed = set()
        else:
            pop = np.array(init, float); fixed = set(range(1, 1 + len(seeds))) | ({0} if args.x0 else set())
        pen = rules_pen(pop)
        bad = [i for i in range(len(pop)) if pen[i] > 0 and i not in fixed]
        n0 = len(pop) - int((pen > 0).sum()); ncand = 0; rounds = 0; spare = []
        while bad and rounds < args.init_admissible:
            cand = lo + qmc.LatinHypercube(d=D, seed=rng).random(args.init_pool * NP) * (hi - lo)
            cp = rules_pen(cand); ncand += len(cand); rounds += 1
            ok = cand[cp == 0]
            spare += list(zip(cp[cp > 0], cand[cp > 0]))
            for x in ok[:len(bad)]:
                i = bad.pop(0); pop[i] = x; pen[i] = 0.0
        if bad:                                     # still short: the lowest-penalty candidates seen
            spare.sort(key=lambda t: t[0])
            for (p, x) in spare[:len(bad)]:
                i = bad.pop(0); pop[i] = x; pen[i] = p
            print(f"WARNING admissible init: {int((pen > 0).sum())} members still inadmissible after {rounds} rounds",
                  flush=True)
        na = int((pen == 0).sum())
        print(f"admissible init: {na} of {len(pop)} admissible (before {n0}; seeds {len(fixed)}, rounds {rounds}, "
              f"acceptance {(na - n0) / max(ncand, 1):.4f} of {ncand} candidates)", flush=True)
        init = pop
    ck = None; phaseA = {}
    if args.resume:
        ck = np.load(args.resume); init = ck["x"]; print(f"resume from {args.resume}: {len(init)} members", flush=True)
        phaseA = dict(phaseA_generations=int(ck["nit"]) if "nit" in ck.files else None)
        ja = args.resume[:-len("_ckpt.npz")] + ".json" if args.resume.endswith("_ckpt.npz") else ""
        if ja and os.path.exists(ja):
            fa = json.load(open(ja))
            phaseA.update(phaseA_nfev=fa.get("nfev"), phaseA_minutes=fa.get("minutes"), phaseA_de_fun=fa.get("de_fun"))
        else:
            phaseA.update(phaseA_nfev=None, phaseA_minutes=None, phaseA_de_fun=None)
        print(f"phase A: {phaseA}", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.save)), exist_ok=True)
    from scipy.optimize import differential_evolution
    tg = [time.time()]

    def checkpoint(intermediate_result):
        ir = intermediate_result
        np.savez(args.save + "_ckpt.npz", x=ir.population, f=ir.population_energies, best=ir.x, nit=ir.nit)
        now = time.time()
        print(f"  ckpt nit {ir.nit}: best {ir.fun:.3f}, objective call {_POP.get('sec', 0):.2f}s "
              f"({_POP.get('n', 0)} admissible), generation {now - tg[0]:.2f}s", flush=True)
        tg[0] = now
    t0 = time.time()
    phaseB = dict(mode="de")
    if args.maxiter > 0 and ck is not None and args.phaseb == "rescore":
        xb, fun, nfev, phaseB = rescore_polish(ck, lo, hi, args)
    elif args.maxiter > 0:
        r = differential_evolution(objective_gpu, bounds, maxiter=args.maxiter, popsize=args.popsize, init=init,
                                   polish=False, seed=args.seed, strategy=args.strategy, disp=True, vectorized=True, updating="deferred",
                                   tol=1e-6, callback=checkpoint)
        xb, fun, nfev = r.x, float(r.fun), int(r.nfev)
    else:
        xb, fun, nfev = seeds[0][2], float("nan"), 0
    a, P = unpack5(xb); Pf = {**fil, **P}
    chi2, tab = table(G5, T5, a, Pf, "l5")
    tab.to_csv(args.save + ".csv", index=False)
    out = dict(model="v2", a=a, pre=P, chi2=chi2, de_fun=fun, nfev=nfev, n_targets=len(tab), minutes=(time.time() - t0) / 60,
               args=vars(args), bap_gate=batch_v2.BAP_GATE, objective="gpu_v6_rho", n_free=len(bounds),
               v5=dict(v5_mode=int(Pf.get("v5_mode", 0)), theta_V=float(Pf["theta_V"]), A_eCB=float(Pf.get("A_eCB", 0.0)),
                       theta_eCB=float(Pf.get("theta_eCB", 0.0)), tau_d_NMDA=float(Pf.get("tau_d_NMDA", MV.DEFAULTS["tau_d_NMDA"])),
                       dpre_min=float(Pf["dpre_min"]), tau_E1=float(Pf["tau_E1"]), rho_gamma=rho_v4.opts(Pf)[0]),
               gamma_d=rho_v4.rates(Pf)[0], gamma_p=rho_v4.rates(Pf)[1], fit_gamma=int(bool(args.fit_gamma)),
               phaseB=phaseB, **phaseA)
    print(tab.round(3).to_string(index=False))
    print(f"chi2 L5 {chi2:.2f} over {len(tab)} targets (gamma_d {out['gamma_d']:.3f}, gamma_p {out['gamma_p']:.3f}, "
          f"v5 {out['v5']})")
    tot = chi2
    for G, name in _M[1:]:
        cm, tm = table(G, TM[name], a, Pf, name); tm.to_csv(f"{args.save}_{name}.csv", index=False)
        out[f"chi2_{name}"] = cm; out[f"n_targets_{name}"] = len(tm); tot += cm
        print(tm.round(3).to_string(index=False)); print(f"chi2 {name} {cm:.2f} over {len(tm)} targets")
    out["chi2_total"] = tot
    out["models"] = [dict(name="l5", groups=args.groups, dirs=args.dirs, basis=L5_BASIS)] + [
        dict(name=n, groups=g, dirs=d, basis=b, geom=ge, pairs=p) for n, g, d, b, ge, p in specs]
    k = len(bounds); n = sum(out.get(f"n_targets_{m}", len(tab)) if m != "l5" else len(tab) for _, m in _M)
    out["aic"] = tot + 2 * k; out["bic"] = tot + k * np.log(n)
    if HINGE:        # v7: hinge part of the DE objective at the final point (z from the unweighted csv tables)
        hp = 0.0
        for t, m in [(tab, "l5")] + [(pd_read(f"{args.save}_{m}.csv"), m) for _, m in _M[1:]]:
            for tt, cc, z in zip(t.target, t.condition, t.z):
                for key in (f"{tt}|{cc}", f"{m}/{tt}|{cc}"):
                    if key in HINGE:
                        kk, lam = HINGE[key]
                        hp += float(lam) * (max(abs(z) - float(kk), 0.0) ** 2 if np.isfinite(z) else 1e3); break
        out["hinge"] = {k: list(v) for k, v in HINGE.items()}; out["hinge_pen"] = hp
        print(f"hinge penalty at the final point {hp:.4f} ({out['hinge']})", flush=True)
    if WEIGHTS:      # unweighted chi2 above; split LTD / LTP and weighted objective value of the final point
        tabs = [(tab, "l5")] + [(pd_read(f"{args.save}_{m}.csv"), m) for _, m in _M[1:]]
        z2 = lambda t, msk: float((t.z[msk] ** 2).sum())
        out["weights"] = WEIGHTS
        out["chi2_ltd"] = sum(z2(t, t.ltd) for t, _ in tabs); out["chi2_ltp"] = sum(z2(t, ~t.ltd) for t, _ in tabs)
        out["n_ltd"] = int(sum(t.ltd.sum() for t, _ in tabs)); out["n_ltp"] = int(sum((~t.ltd).sum() for t, _ in tabs))
        out["chi2_weighted"] = float(sum((t.w * t.z ** 2).sum() for t, _ in tabs))
        print(f"UNWEIGHTED chi2 total {tot:.6f} = LTD {out['chi2_ltd']:.3f} ({out['n_ltd']}) + LTP {out['chi2_ltp']:.3f} "
              f"({out['n_ltp']}); weighted objective {out['chi2_weighted']:.3f} (de_fun {fun:.3f})", flush=True)
    print(f"chi2 total {tot:.6f} over {len(_M)} models, {n} targets, k {k}: AIC {out['aic']:.2f}, BIC {out['bic']:.2f}")
    spec6 = dict(scaled=bool(sV or sE or tpre), cexp_dir=args.cexp_dir, scale_V=sV, scale_E=sE, theta_pre=tpre)
    out["v6"] = dict(spec6, k_V=float(Pf["theta_V"]) if sV else None, k_E=float(Pf.get("theta_eCB", 0.0)) if sE else None,
                     seed_medians=med, pool_per_uM=POOL_PER_UM)
    if args.maxiter <= 0 and seeds and seeds[0][1].get("objective") in ("gpu_v5_rho", "gpu_v6_rho"):    # repro of a json
        fj = seeds[0][1]; f6 = fj.get("v6", dict(scaled=False, cexp_dir="", scale_V={}, scale_E={}, theta_pre={}))
        same = all(f6.get(k) == v for k, v in spec6.items() if k != "cexp_dir") and (not spec6["scaled"] or
                                                                                   f6.get("cexp_dir") == spec6["cexp_dir"])
        if not same:
            print(f"repro skipped: seed rule spec {f6} differs from this run {spec6}", flush=True)
        else:
            ok = True; dmax = 0.0
            for _, name in _M:
                kk = "chi2" if name == "l5" else f"chi2_{name}"
                if kk in fj:
                    d = abs(out[kk] - fj[kk]); ok &= d <= 1e-6; dmax = max(dmax, d)
                    print(f"repro {name}: json {fj[kk]:.12f}  now {out[kk]:.12f}  diff {d:.1e}", flush=True)
            d = abs(out["chi2_total"] - fj.get("chi2_total", np.nan)); ok &= bool(d <= 1e-6); dmax = max(dmax, d)
            print(f"repro total: json {fj.get('chi2_total', np.nan):.12f}  now {out['chi2_total']:.12f}  diff {d:.1e}",
                  flush=True)
            print(f"REPRO OK (1e-6; max diff {dmax:.1e}{', bit-identical' if dmax == 0 else ''})" if ok else "REPRO DIFF",
                  flush=True)
            if not ok:
                json.dump(out, open(args.save + ".json", "w"), indent=1)
                raise SystemExit("REPRO DIFF")
    json.dump(out, open(args.save + ".json", "w"), indent=1)
    print(f"-> {args.save}.json")


if __name__ == "__main__":
    main()
