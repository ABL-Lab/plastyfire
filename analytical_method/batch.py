"""All synapses in one array, so evaluating an a-param set is a few numpy ops.

Two evaluation modes:

  full  — steps the rho ODE over every one of the 21001 decimated samples.
          Reference accuracy, ~seconds per parameter set.
  fast  — mean-field. tau_ind is 70 s against a 42 s protocol, so rho barely
          moves during any single calcium transient and only the FRACTION of
          time above each threshold matters. That fraction comes from the
          per-synapse sorted effcai by binary search, making a parameter set
          cost a searchsorted plus a short scalar ODE. Must be checked against
          `full` before it is trusted (see validate.py).
"""
import glob, os
import numpy as np, pandas as pd

from constants import (RHO_STAR_GB, TAU_IND_GB, GAMMA_D_GB, GAMMA_P_GB,
                       RHO_BINARY_THRESHOLD, PROTO_DT)
import model as M

HERE = os.path.dirname(os.path.abspath(__file__))
# Override with ANALYTICAL_EXTRACTED to point at a different decimation, e.g.
#   ANALYTICAL_EXTRACTED=extracted_d0p25 python validate.py
EXTRACTED = os.path.join(HERE, os.environ.get("ANALYTICAL_EXTRACTED", "extracted"))


class Batch:
    def __init__(self, limit=None, extracted=None, stride=1):
        """stride: decimate effcai on load. The dt scan (dt_scan.py) gives
        0.25 ms -> maxcurve 0.00792 (the trace-interpolation floor, 10.3 s/eval)
        1.00 ms -> maxcurve 0.00875 at 4x less time and RAM.
        So fit at stride 4 from the 0.25 ms data, validate at stride 2."""
        src = extracted or EXTRACTED
        fs = sorted(glob.glob(os.path.join(src, "*.npz")))
        if not fs:
            raise FileNotFoundError(f"no .npz in {src}")
        if limit: fs = fs[:limit]
        eff, cp, cq, r0, ro, grp, meta = [], [], [], [], [], [], []
        for gi, f in enumerate(fs):
            pair, proto = os.path.basename(f)[:-4].split("__")
            d = np.load(f)
            n = len(d["syn"])
            e = d["effcai"]
            eff.append(e[:, ::stride] if stride > 1 else e)
            cp.append(d["c_pre"]); cq.append(d["c_post"])
            r0.append(d["rho0"]); ro.append(d["rho_obs"]); grp.append(np.full(n, gi))
            meta.append(dict(pair=pair, proto=proto, dt=PROTO_DT[proto],
                             n=n, dt_ms=float(d["dt_ms"]) * stride))
        self.effcai = np.concatenate(eff)                  # (N, T) float32
        # (T, N) C-contiguous copy: rho_full touches one timestep at a time, and
        # a column slice of (N,T) is strided. This is the inner loop of every
        # optimiser evaluation, so the extra copy pays for itself immediately.
        self.effcaiT = np.ascontiguousarray(self.effcai.T)
        self.c_pre  = np.concatenate(cp);  self.c_post = np.concatenate(cq)
        self.rho0   = np.concatenate(r0);  self.rho_obs = np.concatenate(ro)
        self.group  = np.concatenate(grp)
        self.meta   = meta
        self.dt_ms  = meta[0]["dt_ms"]
        self.T      = self.effcai.shape[1]
        self.N      = self.effcai.shape[0]
        self.peak   = self.effcai.max(axis=1).astype(np.float64)  # per-synapse max effcai
        self.src    = src
        self._sorted = None      # built lazily; doubles memory when it is

    @property
    def sorted(self):
        if self._sorted is None:
            self._sorted = np.sort(self.effcai, axis=1)
        return self._sorted

    # ---- theta -----------------------------------------------------------
    def thetas(self, a, is_apical=None):
        ap = np.zeros(self.N, bool) if is_apical is None else is_apical
        d0 = np.where(ap, a["a20"], a["a00"]); d1 = np.where(ap, a["a21"], a["a01"])
        p0 = np.where(ap, a["a30"], a["a10"]); p1 = np.where(ap, a["a31"], a["a11"])
        return (d0 * self.c_pre + d1 * self.c_post,
                p0 * self.c_pre + p1 * self.c_post)

    # ---- rho -------------------------------------------------------------
    def rho_full(self, td, tp, gamma_d=GAMMA_D_GB, gamma_p=GAMMA_P_GB):
        h = (self.dt_ms / 1000.0) / TAU_IND_GB
        rho = self.rho0.astype(np.float64).copy()
        E = self.effcaiT
        for k in range(self.T):
            e = E[k]
            dep = e > td; pot = e > tp
            # (1 - pot) : depression is GATED OFF while potentiation is active.
            # This factor is in the compiled mechanism the sims load
            # (DEES_cell_packages/other_mods/modified_mechanisms/GluSynapse.mod)
            # but NOT in plastyfire/GluSynapse.mod. Omitting it made every
            # synapse depress above theta_p and biased rho low by ~0.13.
            rho += h * (-rho * (1 - rho) * (RHO_STAR_GB - rho)
                        + pot * gamma_p * (1 - rho)
                        - dep * (1 - pot) * gamma_d * rho)
            np.clip(rho, 0.0, 1.0, out=rho)
        return rho

    def _frac_above(self, th):
        """Fraction of samples with effcai > th, per synapse, via binary search."""
        S = self.sorted
        idx = np.array([np.searchsorted(S[i], th[i], side="right")
                        for i in range(self.N)])
        return 1.0 - idx / self.T

    def rho_fast(self, td, tp, gamma_d=GAMMA_D_GB, gamma_p=GAMMA_P_GB, steps=400):
        fd = self._frac_above(td); fp = self._frac_above(tp)
        h = (self.dt_ms * self.T / 1000.0) / TAU_IND_GB / steps
        rho = self.rho0.astype(np.float64).copy()
        for _ in range(steps):
            # mean-field analogue of dep*(1-pot): depression only acts during the
            # time in the band (theta_d, theta_p), i.e. fd - fp, not fd.
            rho += h * (-rho * (1 - rho) * (RHO_STAR_GB - rho)
                        + fp * gamma_p * (1 - rho)
                        - np.maximum(fd - fp, 0.0) * gamma_d * rho)
            np.clip(rho, 0.0, 1.0, out=rho)
        return rho

    # ---- EPSP ------------------------------------------------------------
    def curve(self, rho_final):
        rows = []
        off = 0
        for gi, m in enumerate(self.meta):
            sl = slice(off, off + m["n"]); off += m["n"]
            pre, post = m["pair"].split("-")
            try:
                r = M.epsp_ratio(M.load_basis(pre, post),
                                 self.rho0[sl], rho_final[sl])
            except Exception:
                continue
            if r is None: continue
            rows.append(dict(pair=m["pair"], dt=m["dt"], ratio=r))
        df = pd.DataFrame(rows)
        summ = (df.groupby("dt")["ratio"].agg(["mean", "sem", "count"])
                  .reset_index().sort_values("dt"))
        return summ, df

    def evaluate(self, a, mode="fast", is_apical=None,
                 gamma_d=GAMMA_D_GB, gamma_p=GAMMA_P_GB):
        td, tp = self.thetas(a, is_apical)
        f = self.rho_fast if mode == "fast" else self.rho_full
        rho = f(td, tp, gamma_d, gamma_p)
        return rho, self.curve(rho)
