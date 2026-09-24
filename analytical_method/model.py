"""Stage 2 (cheap): a-params -> theta -> rho -> EPSP ratio -> STDP curve.

    theta_d = a00*c_pre + a01*c_post      (basal)    a20/a21 (apical)
    theta_p = a10*c_pre + a11*c_post      (basal)    a30/a31 (apical)

    rho' = [ -rho(1-rho)(rho* - rho)
             + pot*gamma_p*(1-rho) - dep*gamma_d*rho ] / (1e3*tau_ind)

with dep = 1 while effcai > theta_d, pot = 1 while effcai > theta_p. tau_ind is
70 s against a 42 s protocol, so rho moves slowly and an explicit step on the
decimated grid is stable.

EPSP uses the basis exactly as evaluator_edges does — linear superposition over
synapses whose FINAL rho >= 0.5. That binarisation is why the objective is
piecewise-constant in the a-params, and why a gradient method cannot be used.
"""
import os, sys, glob
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from constants import (RHO_STAR_GB, TAU_IND_GB, GAMMA_D_GB, GAMMA_P_GB,
                       RHO_BINARY_THRESHOLD, PROTO_DT)

ROOT      = "/lustre06/project/6077694/dhuruva/plastyfire"
BASIS_DIR = os.path.join(ROOT, "basis_results_edges_ion_channels")
EXTRACTED = os.path.join(os.path.dirname(os.path.abspath(__file__)), "extracted")


def thetas(c_pre, c_post, a, is_apical):
    """Per-synapse (theta_d, theta_p). is_apical selects a2x/a3x over a0x/a1x."""
    d0, d1 = np.where(is_apical, a["a20"], a["a00"]), np.where(is_apical, a["a21"], a["a01"])
    p0, p1 = np.where(is_apical, a["a30"], a["a10"]), np.where(is_apical, a["a31"], a["a11"])
    return d0 * c_pre + d1 * c_post, p0 * c_pre + p1 * c_post


def integrate_rho(effcai, dt_ms, theta_d, theta_p, rho0):
    """Vectorised over synapses. effcai (n_syn, n_t); returns final rho (n_syn,)."""
    h   = (dt_ms / 1000.0) / TAU_IND_GB
    rho = np.asarray(rho0, dtype=np.float64).copy()
    td  = theta_d[:, None]; tp = theta_p[:, None]
    dep = (effcai > td)
    pot = (effcai > tp)
    for k in range(effcai.shape[1]):
        rho += h * (-rho * (1 - rho) * (RHO_STAR_GB - rho)
                    + pot[:, k] * GAMMA_P_GB * (1 - rho)
                    - dep[:, k] * (1 - pot[:, k]) * GAMMA_D_GB * rho)
        np.clip(rho, 0.0, 1.0, out=rho)
    return rho


_BASIS = {}
def load_basis(pre_gid, post_gid):
    key = (pre_gid, post_gid)
    if key not in _BASIS:
        _BASIS[key] = pd.read_csv(os.path.join(BASIS_DIR, f"basis_{pre_gid}_{post_gid}.csv"))
    return _BASIS[key]


def epsp_from_basis(basis_df, rho_vec):
    """Linear superposition with binary rho — mirrors evaluator_edges."""
    n   = len(rho_vec)
    rb  = [1 if r >= RHO_BINARY_THRESHOLD else 0 for r in rho_vec]
    z   = ",".join(["0"] * n)
    r0  = basis_df.loc[basis_df["config"] == z]
    if r0.empty:
        raise ValueError(f"basis missing all-zeros config for n={n}")
    e0, var0 = float(r0["mean"].values[0]), float(r0["std"].values[0]) ** 2
    mean, var, k = e0, 0.0, 0
    for i, b in enumerate(rb):
        if not b: continue
        s = ["0"] * n; s[i] = "1"
        row = basis_df.loc[basis_df["config"] == ",".join(s)]
        if row.empty:
            raise ValueError(f"basis missing singleton {i} for n={n}")
        mean += float(row["mean"].values[0]) - e0
        var  += float(row["std"].values[0]) ** 2
        k    += 1
    var += (1 - k) ** 2 * var0
    return mean, float(np.sqrt(max(0.0, var)))


def epsp_ratio(basis_df, rho_initial, rho_final):
    """Same ratio definition as plot_stdp_ic_best.py, including the cv2 term."""
    b_m, b_s = epsp_from_basis(basis_df, rho_initial)
    a_m, _   = epsp_from_basis(basis_df, rho_final)
    if b_m == 0: return None
    cv2 = min((b_s / b_m) ** 2, 0.25)
    return (a_m / b_m) * (1.0 + cv2)


def load_extracted(limit_pairs=None):
    out = []
    for f in sorted(glob.glob(os.path.join(EXTRACTED, "*.npz"))):
        pair, proto = os.path.basename(f)[:-4].split("__")
        if limit_pairs and pair not in limit_pairs: continue
        d = np.load(f)
        out.append(dict(pair=pair, proto=proto, dt=PROTO_DT[proto],
                        syn=d["syn"], effcai=d["effcai"], dt_ms=float(d["dt_ms"]),
                        rho0=d["rho0"], rho_obs=d["rho_obs"],
                        c_pre=d["c_pre"], c_post=d["c_post"]))
    return out


def stdp_curve(datasets, a, apical_map=None):
    """-> DataFrame[dt, mean, sem, n] plus the per-record detail."""
    rows = []
    for r in datasets:
        is_ap = (np.zeros(len(r["syn"]), bool) if apical_map is None
                 else np.array([apical_map.get(int(s), False) for s in r["syn"]]))
        td, tp = thetas(r["c_pre"], r["c_post"], a, is_ap)
        rho_f  = integrate_rho(r["effcai"], r["dt_ms"], td, tp, r["rho0"])
        pre, post = r["pair"].split("-")
        try:
            ratio = epsp_ratio(load_basis(pre, post), r["rho0"], rho_f)
        except Exception:
            continue
        if ratio is None: continue
        rows.append(dict(pair=r["pair"], dt=r["dt"], ratio=ratio,
                         rho_pred=rho_f.mean(), rho_obs=r["rho_obs"].mean()))
    df = pd.DataFrame(rows)
    if df.empty: return df, df
    summ = (df.groupby("dt")["ratio"].agg(["mean", "sem", "count"])
              .reset_index().sort_values("dt"))
    return summ, df
